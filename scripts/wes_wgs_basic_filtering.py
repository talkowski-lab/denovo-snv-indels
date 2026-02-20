###
# Adapted from Kyle Satterstrom's WES and Shan Dong's WGS filtering scripts.
###

import argparse
import datetime
import hail as hl
import numpy as np
import os
import pandas as pd


# Define sex-aware variant call rate calculation with Hardy-Weinberg p value
def sex_aware_variant_annotations_with_pHWE(mt):
    num_males = mt.aggregate_cols(hl.agg.count_where(mt.is_female == False))
    num_females = mt.aggregate_cols(hl.agg.count_where(mt.is_female == True))

    mt = mt.annotate_rows(
        male_hets=hl.agg.count_where(mt.GT.is_het() & (mt.is_female == False)),
        male_homvars=hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == False)),
        male_calls=hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == False)),
        female_hets=hl.agg.count_where(mt.GT.is_het() & (mt.is_female == True)),
        female_homvars=hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == True)),
        female_calls=hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == True)),
    )

    mt = mt.annotate_rows(
        call_rate=(
            hl.case()
            .when(mt.locus.in_y_nonpar(), (mt.male_calls / num_males))
            .when(
                mt.locus.in_x_nonpar(),
                (mt.male_calls + 2 * mt.female_calls) / (num_males + 2 * num_females),
            )
            .default((mt.male_calls + mt.female_calls) / (num_males + num_females))
        ),
        sex_aware_AC=(
            hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_homvars)
            .when(
                mt.locus.in_x_nonpar(),
                mt.male_homvars + mt.female_hets + 2 * mt.female_homvars,
            )
            .default(
                mt.male_hets
                + 2 * mt.male_homvars
                + mt.female_hets
                + 2 * mt.female_homvars
            )
        ),
        sex_aware_AN=(
            hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_calls)
            .when(mt.locus.in_x_nonpar(), mt.male_calls + 2 * mt.female_calls)
            .default(2 * mt.male_calls + 2 * mt.female_calls)
        ),
        pHWE=(
            hl.case()
            .when(mt.locus.in_y_nonpar() | mt.locus.in_mito(), 1.0)
            .when(
                mt.locus.in_x_nonpar(),
                hl.hardy_weinberg_test(
                    hl.int32(mt.female_calls - mt.female_hets - mt.female_homvars),
                    hl.int32(mt.female_hets),
                    hl.int32(mt.female_homvars),
                ).p_value,
            )
            .default(
                hl.hardy_weinberg_test(
                    hl.int32(
                        mt.male_calls
                        + mt.female_calls
                        - mt.male_hets
                        - mt.female_hets
                        - mt.male_homvars
                        - mt.female_homvars
                    ),
                    hl.int32(mt.male_hets + mt.female_hets),
                    hl.int32(mt.male_homvars + mt.female_homvars),
                ).p_value
            )
        ),
    )
    return mt


# Define sex-aware sample call rate calculation
def sex_aware_sample_annotations(mt):
    num_y_non_par_vars = mt.aggregate_rows(hl.agg.count_where(mt.locus.in_y_nonpar()))
    num_all_other_vars = mt.aggregate_rows(hl.agg.count_where(~mt.locus.in_y_nonpar()))

    mt = mt.annotate_cols(
        sample_call_rate=(
            hl.case()
            .when(
                mt.is_female == True,
                hl.agg.count_where(hl.is_defined(mt.GT) & ~mt.locus.in_y_nonpar())
                / num_all_other_vars,
            )
            .default(
                hl.agg.count_where(hl.is_defined(mt.GT))
                / (num_y_non_par_vars + num_all_other_vars)
            )
        )
    )

    return mt


def main(args):
    hl.init(
        min_block_size=128,
        spark_conf={
            "spark.executor.cores": args.cores,
            "spark.executor.memory": f"{int(np.floor(args.mem*0.4))}g",
            "spark.driver.cores": args.cores,
            "spark.driver.memory": f"{int(np.floor(args.mem*0.4))}g",
        },
        tmp_dir="tmp",
        local_tmpdir="tmp",
    )

    bucket_id = args.bucket_id
    prefix = os.path.basename(args.annot_mt).split("_wes_denovo_annot.mt")[0]

    # Read in MT
    mt = hl.read_matrix_table(args.annot_mt)

    # WES & WGS: Filter to samples in both PED and MT
    tmp_ped = pd.read_csv(args.ped_uri, sep="\t")
    # check tmp_ped number of columns
    if len(tmp_ped.columns) > 6:
        tmp_ped = tmp_ped.iloc[:, :6]
    # Get samples in both PED and MT
    samps = mt.s.collect_as_set().intersection(tmp_ped.iloc[:, 1])
    # Subset MT to these samples
    mt = mt.filter_cols(hl.literal(samps).contains(mt.s))
    # Subset PED to these samples
    tmp_ped = tmp_ped[tmp_ped.iloc[:, 1].isin(samps)]  # sample_id
    tmp_ped = tmp_ped.drop_duplicates(tmp_ped.columns[1])
    tmp_ped.to_csv(f"{prefix}.ped", sep="\t", index=False)
    ped_uri_processed = f"{prefix}.ped"
    ped = hl.import_table(ped_uri_processed, impute=True, delimiter="\t")
    original_cols = list(ped.row.keys())
    new_cols = [
        "family_id",
        "sample_id",
        "paternal_id",
        "maternal_id",
        "sex",
        "phenotype",
    ]
    ped = ped.rename({old: new for old, new in zip(original_cols, new_cols)})
    ped = ped.key_by("sample_id")

    # WES & WGS: Filter low complexity regions
    try:
        lcr = hl.import_bed(args.lcr_uri, reference_genome="GRCh38")
    except:
        lcr = hl.import_bed(args.lcr_uri, reference_genome="GRCh38", force_bgz=True)
    if args.genome_build == "GRCh37":
        rg37 = hl.get_reference("GRCh37")
        rg38 = hl.get_reference("GRCh38")
        rg38.add_liftover(
            "gs://hail-common/references/grch38_to_grch37.over.chain.gz", rg37
        )
        lcr = lcr.annotate(new_locus=hl.liftover(lcr.interval, args.genome_build))
        lcr = lcr.filter(hl.is_defined(lcr.new_locus))
        lcr = lcr.key_by(locus=lcr.new_locus)
    mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)

    # WES & WGS: PASS filters
    if args.filter_snv_pass:
        mt = mt.filter_rows(
            hl.is_snp(mt.alleles[0], mt.alleles[1])
            & ~((hl.len(mt.filters) == 0) | hl.is_missing(mt.filters)),
            keep=False,
        )
    if args.filter_indel_pass:
        mt = mt.filter_rows(
            hl.is_indel(mt.alleles[0], mt.alleles[1])
            & ~((hl.len(mt.filters) == 0) | hl.is_missing(mt.filters)),
            keep=False,
        )

    mt = mt.annotate_entries(DPC=hl.sum(mt.AD), AB=mt.AD[1] / hl.sum(mt.AD))

    # WGS (should also be applied to WES, co-opts WES DP filters): DPC filters
    mt = mt.filter_entries(mt.DPC < args.min_dpc, keep=False)
    mt = mt.filter_entries(mt.DPC > args.max_dpc, keep=False)

    # TODO: Figure out how to adjust for XO
    # WES (should also be applied to WGS): Sex-specific genotype filtering:
    # - Any call in a female with depth less than minimum threshold
    # - Any call on the Y in females
    # - Any autosomal or PAR call in a male with depth less than minimum threshold
    # - Het calls in males in hemizygous regions
    mt = mt.annotate_cols(reported_sex=ped[mt.s].sex)
    mt = mt.annotate_cols(is_female=(mt.reported_sex == 2))
    mt = mt.filter_entries(
        ((mt.is_female == True) & (mt.DP < args.female_min_dp))
        | ((mt.is_female == True) & (mt.locus.in_y_nonpar()))
        | (
            (mt.is_female == False)
            & (mt.locus.in_autosome_or_par())
            & (mt.DP < args.male_auto_min_dp)
        )
        | (
            (mt.is_female == False)
            & (mt.GT.is_het())
            & (mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar())
        ),
        keep=False,
    )

    # WGS: Apply row/variant-level filters
    if args.qual_threshold is not None:
        mt = mt.filter_rows(mt.qual >= args.qual_threshold)
    if args.mq_threshold is not None:
        mt = mt.filter_rows(mt.info.MQ >= args.mq_threshold)

    # WGS: Set SNV row-level variant INFO filters
    snv_cond_row = hl.is_snp(mt.alleles[0], mt.alleles[1])
    if args.sor_threshold_snv is not None:
        snv_cond_row = snv_cond_row & (mt.info.SOR <= args.sor_threshold_snv)
    if args.readposranksum_threshold_snv is not None:
        snv_cond_row = snv_cond_row & (
            mt.info.ReadPosRankSum >= args.readposranksum_threshold_snv
        )
    if args.qd_threshold_snv is not None:
        snv_cond_row = snv_cond_row & (mt.info.QD >= args.qd_threshold_snv)
    # WGS: Set indel row-level variant INFO filters
    indel_cond_row = hl.is_indel(mt.alleles[0], mt.alleles[1])
    if args.sor_threshold_indel is not None:
        indel_cond_row = indel_cond_row & (mt.info.SOR <= args.sor_threshold_indel)
    if args.readposranksum_threshold_indel is not None:
        indel_cond_row = indel_cond_row & (
            mt.info.ReadPosRankSum >= args.readposranksum_threshold_indel
        )
    if args.qd_threshold_indel is not None:
        indel_cond_row = indel_cond_row & (mt.info.QD >= args.qd_threshold_indel)
    # WGS: Apply row-level variant INFO filters
    mt = mt.filter_rows(snv_cond_row | indel_cond_row, keep=True)

    # WGS (should also be applied to WES, co-opts a WES Het AB filter): HomRef, Het AB filters
    homref_snv_indel_cond = mt.GT.is_hom_ref() & (mt.AB > args.max_homref_ab)
    het_snv_cond = (
        hl.is_snp(mt.alleles[0], mt.alleles[1])
        & mt.GT.is_het()
        & ((mt.AB < args.min_het_snv_ab) | (mt.AB > args.max_het_snv_ab))
    )
    het_indel_cond = (
        hl.is_indel(mt.alleles[0], mt.alleles[1])
        & mt.GT.is_het()
        & ((mt.AB < args.min_het_indel_ab) | (mt.AB > args.max_het_indel_ab))
    )
    mt = mt.filter_entries(
        homref_snv_indel_cond | het_snv_cond | het_indel_cond, keep=False
    )

    # WES (should also be applied to WES): Filter by het pAB of 0
    mt = mt.filter_entries(
        mt.GT.is_het() & (mt.pAB < args.het_pab_threshold), keep=False
    )

    # WGS: Variant GQ mean filter
    if args.min_mean_gq is not None and not args.exclude_gq_filters:
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= args.min_mean_gq, keep=True)

    # WGS (should also be applied to WES, co-opts a WES HomRef GQ filter): GQ filters
    if args.min_gq is not None and not args.exclude_gq_filters:
        mt = mt.filter_entries(mt.GQ < args.min_gq, keep=False)

    # WES: HomVar and Het PL filter
    if args.min_pl is not None:
        mt = mt.filter_entries(
            (mt.GT.is_hom_var() | mt.GT.is_het()) & (mt.PL[0] < args.min_pl),
            keep=False,
        )

    # WES: Filter by informative read threshold
    if args.informative_read_threshold is not None:
        mt = mt.filter_entries(
            mt.GT.is_het() & (mt.DPC < (args.informative_read_threshold * mt.DP)),
            keep=False,
        )
        mt = mt.filter_entries(
            mt.GT.is_hom_var() & (mt.AD[1] < (args.informative_read_threshold * mt.DP)),
            keep=False,
        )

    # Add sex-aware variant annotations including pHWE
    mt = sex_aware_variant_annotations_with_pHWE(mt)

    # WES: Variant call rate filter
    if args.call_rate_threshold is not None:
        mt = mt.filter_rows(mt.call_rate < args.call_rate_threshold, keep=False)
    # WES (should also apply to WGS): pHWE filter
    mt = mt.filter_rows(mt.pHWE < args.phwe_threshold, keep=False)

    # Drop variants with no alleles in sample set
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep=True)
    mt = mt.drop("variant_qc")

    # Calculate sample call rates
    mt = sex_aware_sample_annotations(mt)

    # get sample-level stats to plot
    hl.sample_qc(mt).cols().flatten().export(
        f"{prefix}_final_annot_post_filter_qc_info.txt"
    )

    # export mt
    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{prefix}_basic_filtering.mt"
    pd.Series([filename]).to_csv("mt_uri.txt", index=False, header=None)
    mt.write(filename, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Basic filtering for WES/WGS rare/de novo calling"
    )
    parser.add_argument("--annot_mt", required=True, help="Input MatrixTable")
    parser.add_argument("--ped_uri", required=True)
    parser.add_argument("--cores", default="8")
    parser.add_argument("--mem", type=float, required=True, help="Memory in GB")
    parser.add_argument("--bucket_id", required=True)
    parser.add_argument("--lcr_uri", required=True)
    parser.add_argument("--genome_build", default="GRCh38")

    # Parameters for hardcoded filters
    parser.add_argument("--filter_snv_pass", action="store_true")
    parser.add_argument("--filter_indel_pass", action="store_true")
    parser.add_argument("--exclude_gq_filters", action="store_true")
    parser.add_argument("--min_dpc", type=int, required=True)
    parser.add_argument("--max_dpc", type=int, required=True)
    parser.add_argument("--female_min_dp", type=int, required=True)
    parser.add_argument("--male_auto_min_dp", type=int, required=True)
    parser.add_argument("--qual_threshold", type=int)
    parser.add_argument("--mq_threshold", type=int)
    parser.add_argument("--sor_threshold_snv", type=float)
    parser.add_argument("--readposranksum_threshold_snv", type=float)
    parser.add_argument("--qd_threshold_snv", type=float)
    parser.add_argument("--sor_threshold_indel", type=float)
    parser.add_argument("--readposranksum_threshold_indel", type=float)
    parser.add_argument("--qd_threshold_indel", type=float)
    parser.add_argument("--max_homref_ab", type=float, required=True)
    parser.add_argument("--min_het_snv_ab", type=float, required=True)
    parser.add_argument("--max_het_snv_ab", type=float, required=True)
    parser.add_argument("--min_het_indel_ab", type=float, required=True)
    parser.add_argument("--max_het_indel_ab", type=float, required=True)
    parser.add_argument("--max_het_indel_ab", type=float, required=True)
    parser.add_argument("--het_pab_threshold", type=float, required=True)
    parser.add_argument("--min_mean_gq", type=float)
    parser.add_argument("--min_gq", type=int)
    parser.add_argument("--min_pl", type=int)
    parser.add_argument("--informative_read_threshold", type=float)
    parser.add_argument("--call_rate_threshold", type=float)
    parser.add_argument("--phwe_threshold", type=float, required=True)

    args = parser.parse_args()
    main(args)
