import hail as hl

import seaborn as sns
import matplotlib.pyplot as plt

from typing import Union
from tenacity import retry, after_log, before_sleep_log, retry_if_exception_type, stop_after_attempt, wait_exponential
from firecloud import fiss
from firecloud.errors import FireCloudServerError
import firecloud.api as fapi
import re
import sys
import time
import math
import ast

import pandas as pd
import numpy as np
import os

import logging
from logging import INFO, DEBUG
logger = logging.getLogger()
logger.setLevel(INFO)

pd.set_option('display.float_format', '{:.10f}'.format)

import matplotlib.pyplot as plt
import seaborn as sns

import argparse

def parse_args():
    p = argparse.ArgumentParser()

    add = p.add_argument
    add("--vep-vcf-uri", required=True)
    add("--filt-mt-uri", required=True)
    add("--ped-uri", required=True)

    add("--gnomad-non-neuro-af-threshold", type=float, default=0.001)
    add("--cohort-ac-threshold", type=int, default=20)
    add("--cohort-af-threshold", type=float, default=0.001)
    add("--affected-ac-threshold", type=int, default=None)  # Optional
    add("--affected-af-threshold", type=float, default=None)  # Optional
    add("--coding-only", action="store_true")
    add("--mem", type=float, default=4)

    return p.parse_args()


args = parse_args()

vep_vcf_uri = args.vep_vcf_uri
filt_mt_uri = args.filt_mt_uri
ped_uri = args.ped_uri
gnomad_non_neuro_af_threshold = args.gnomad_non_neuro_af_threshold
cohort_ac_threshold = args.cohort_ac_threshold
cohort_af_threshold = args.cohort_af_threshold
affected_ac_threshold = args.affected_ac_threshold
affected_af_threshold = args.affected_af_threshold
coding_only = args.coding_only
mem = args.mem

prefix = os.path.basename(filt_mt_uri).split('.mt')[0]

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.executor.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'}, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

# Functions

# Parent-aware TDT annotations
#
# This draws heavily from Jack Kosmicki's function, 
# https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test,
# and, when summed, should give the same output values
#
# Requires trio dataset and, like Hail's TDT function, does not cover the Y
#
# Also note:
# 1) Uses sexes from ped and assumes fathers are male and mothers are female
# 2) Does not allow fathers to be het when they should be hemizygous
# 3) To match Jack's function, requires fathers to have a genotype even when considering regions where 
#    proband is hemizygous
def parent_aware_t_u_annotations_v4(td):

    # First decide copy state
    td = td.annotate_entries(autosomal_copy_state = ( hl.case()
            .when(td.locus.in_autosome() | td.locus.in_x_par() | (td.is_female == True), True)
            .when(td.locus.in_x_nonpar() & (td.is_female == False), False)
            .default(hl.missing('bool')) ) )
    # Note: the above uses the "is_female" from the ped and not from the dataset itself
    
    # Now annotate t & u values
    td = td.annotate_entries(
        
        t_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_var()) )) ), 1, 0),
        
        t_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_var()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_var() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),
        # I could consider removing any reference at all to father's genotype in this last line    
        
        u_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_ref()) )) ), 1, 0),
        
        u_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_ref()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_ref() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),   
        # Again, could consider removing any reference at all to father's genotype in this last line
        
        t_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0),

        u_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0)        
    )
        
    return (td)

# try this from https://github.com/Nealelab/recessive/blob/401af812dc4f1fc51cf2e8912aa598e2cef44e3c/Hail_%26_Export_Pipeline_Genotyped_dataset.ipynb
from typing import *

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost"]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant"
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant"]

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant"
]

coding_variants = ['coding_sequence_variant', 'frameshift_variant', 
        'incomplete_terminal_codon_variant', 'inframe_deletion', 'inframe_insertion',
        'missense_variant', 'protein_altering_variant', 'splice_acceptor_variant',
        'splice_donor_variant', 'start_lost', 'stop_gained', 'stop_lost',
        'stop_retained_variant', 'synonymous_variant']

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING + ['.', 'NA']

def filter_vep_to_canonical_transcripts(mt: Union[hl.MatrixTable, hl.Table],
                                        vep_root: str = 'vep') -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(lambda csq: csq.CANONICAL == "YES")
    vep_data = mt[vep_root].annotate(transcript_consequences=hl.if_else(canonical.size()>0, canonical, mt[vep_root].transcript_consequences))
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

def add_most_severe_consequence_to_consequence(tc: hl.expr.StructExpression) -> hl.expr.StructExpression:

    """
    Add most_severe_consequence annotation to transcript consequences
    This is for a given transcript, as there are often multiple annotations for a single transcript:
    e.g. splice_region_variant&intron_variant -> splice_region_variant
    """

    csqs = hl.literal(CSQ_ORDER)

    return tc.annotate(
        most_severe_consequence=csqs.find(lambda c: tc.Consequence.contains(c))
    )

def process_consequences(mt: Union[hl.MatrixTable, hl.Table], vep_root: str = 'vep',
                         penalize_flags: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq_by_gene, any_LoF into [vep_root]
    :param MatrixTable mt: Input MT
    :param str vep_root: Root for vep annotation (probably vep)
    :param bool penalize_flags: Whether to penalize LoFTEE flagged variants, or treat them as equal to HC
    :return: MT with better formatted consequences
    :rtype: MatrixTable
    """
    csqs = hl.literal(CSQ_ORDER)
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def find_worst_transcript_consequence(tcl: hl.expr.ArrayExpression) -> hl.expr.StructExpression:
        """
        Gets worst transcript_consequence from an array of em
        """
        flag_score = 500
        no_flag_score = flag_score * (1 + penalize_flags)
        non_coding_score = 600
        non_canonical_score = 500
        def csq_score(tc):
#             return csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
            return csq_dict.get(tc.most_severe_consequence, csq_dict['NA'])

        # EDITED: PENALIZE NON-CODING AND NON-CANONICAL
        tcl = tcl.map(lambda tc: tc.annotate(
            csq_score=hl.case(missing_false=True)
            .when((tc.BIOTYPE != 'protein_coding'), csq_score(tc) + non_coding_score)
            .when((tc.CANONICAL != 'YES'), csq_score(tc) + non_canonical_score)
            # .when((tc.LoF == 'HC') & (tc.LoF_flags == ''), csq_score(tc) - no_flag_score)
            # .when((tc.LoF == 'HC') & (tc.LoF_flags != ''), csq_score(tc) - flag_score)
            # .when(tc.LoF == 'LC', csq_score(tc) - 10) 
            .when(tc.PolyPhen.contains('probably_damaging'), csq_score(tc) - 0.5)  # EDITED
            .when(tc.PolyPhen.contains('possibly_damaging'), csq_score(tc) - 0.25)
            .when(tc.PolyPhen.contains('benign'), csq_score(tc) - 0.1)
            .default(csq_score(tc))
        ))
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    transcript_csqs = mt[vep_root].transcript_consequences.map(add_most_severe_consequence_to_consequence)

    gene_dict = transcript_csqs.group_by(lambda tc: tc.SYMBOL)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)
    lowest_score = hl.or_missing(hl.len(sorted_scores) > 0, sorted_scores[0].csq_score)
    gene_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.SYMBOL)
    ensg_with_worst_csq = sorted_scores.filter(lambda tc: tc.csq_score == lowest_score).map(lambda tc: tc.Gene)

    vep_data = mt[vep_root].annotate(transcript_consequences=transcript_csqs,
                                     worst_consequence_term=csqs.find(lambda c: transcript_csqs.map(lambda csq: csq.most_severe_consequence).contains(c)),
                                     worst_csq_by_gene=sorted_scores,  # EDITED
                                     worst_csq=sorted_scores[0],
                                    #  any_LoF=hl.any(lambda x: x.LoF == 'HC', worst_csq_gene),
                                     gene_with_most_severe_csq=gene_with_worst_csq,
                                     ensg_with_most_severe_csq=ensg_with_worst_csq)
    
    return mt.annotate_rows(**{vep_root: vep_data}) if isinstance(mt, hl.MatrixTable) else mt.annotate(**{vep_root: vep_data})

def apply_cohort_ac_af_filters(mt, cohort_ac_threshold, cohort_af_threshold):
    return mt.filter_rows(
        (mt.info.cohort_AC <= cohort_ac_threshold) |
        (mt.info.cohort_AF <= cohort_af_threshold)
    )

def apply_affected_ac_filter(mt, affected_ac_threshold):
    return mt.filter_rows(mt.affected_AC <= affected_ac_threshold)

def annotate_affected_unaffected_AC(mt, ped_ht):
    # Annotate phenotype in MT
    mt = mt.annotate_cols(phenotype=ped_ht[mt.s].phenotype)
    
    # Pre-calculate counts of individuals to determine AN
    n_unaffected = ped_ht.filter(ped_ht.phenotype == 1).count()
    n_affected = ped_ht.filter(ped_ht.phenotype == 2).count()

    # Get cohort unaffected/affected het and homvar counts
    mt = mt.annotate_rows(**{
        "n_het_unaffected": hl.agg.filter(mt.phenotype==1, hl.agg.sum(mt.GT.is_het())),
        "n_hom_var_unaffected": hl.agg.filter(mt.phenotype==1, hl.agg.sum(mt.GT.is_hom_var())),
        "n_het_affected": hl.agg.filter(mt.phenotype==2, hl.agg.sum(mt.GT.is_het())),
        "n_hom_var_affected": hl.agg.filter(mt.phenotype==2, hl.agg.sum(mt.GT.is_hom_var()))
    })
    
    mt = mt.annotate_rows(
        unaffected_AC = mt.n_het_unaffected + 2*mt.n_hom_var_unaffected,
        affected_AC = mt.n_het_affected + 2*mt.n_hom_var_affected
    )

    # Calculate AF, handling division by zero if a group is empty
    mt = mt.annotate_rows(
        unaffected_AF = hl.if_else(n_unaffected > 0, mt.unaffected_AC / (2 * n_unaffected), 0.0),
        affected_AF = hl.if_else(n_affected > 0, mt.affected_AC / (2 * n_affected), 0.0)
    )

    return mt

# Load
header = hl.get_vcf_metadata(vep_vcf_uri)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

filt_mt = hl.read_matrix_table(filt_mt_uri)

transcript_consequences = filt_mt.info.CSQ.map(lambda csq_str: csq_str.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: '.' if col!='Consequence' else hl.array(['.'])  
                                                        for i, col in enumerate(csq_columns)})))

filt_mt=filt_mt.annotate_rows(vep=hl.Struct(transcript_consequences = transcript_consequences_strs))
filt_mt = filter_vep_to_canonical_transcripts(filt_mt)
filt_mt = process_consequences(filt_mt)
filt_mt = filt_mt.annotate_rows(worst_csq=filt_mt.vep.worst_csq)

# GnomAD AF filter
base_mt = filt_mt.filter_rows(filt_mt.gnomad_non_neuro_AF > gnomad_non_neuro_af_threshold, keep=False)

ped_ht = hl.import_table(ped_uri, delimiter='\t',
                          types={'phenotype': hl.tfloat, 'sex': hl.tfloat}).key_by('sample_id')
base_mt = annotate_affected_unaffected_AC(base_mt, ped_ht)

# Apply affected AC or AF filters if defined
conditions = []
if affected_ac_threshold is not None:
    conditions.append(base_mt.affected_AC <= affected_ac_threshold)
if affected_af_threshold is not None:
    conditions.append(base_mt.affected_AF <= affected_af_threshold)

if conditions:
    case_filtered_mt = base_mt.filter_rows(hl.any(lambda x: x, conditions))
else:
    # Fallback to standard cohort filters if no affected filters defined
    case_filtered_mt = apply_cohort_ac_af_filters(base_mt, cohort_ac_threshold, cohort_af_threshold)

# Apply cohort AC, AF filters
ultra_rare_mt = apply_cohort_ac_af_filters(
    base_mt,
    cohort_ac_threshold,
    cohort_af_threshold
)

ultra_rare_filt_mt_uri = f"{prefix}.ultra.rare.mt"
ultra_rare_mt = ultra_rare_mt.checkpoint(ultra_rare_filt_mt_uri, overwrite=True)

# Pedigree
ped_df = pd.read_csv(ped_uri, sep='\t')
cropped_ped_uri = f"{ped_uri.split('.ped')[0]}_cropped.ped"
ped_df.iloc[:, :6].to_csv(cropped_ped_uri, index=False, sep='\t')

pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')
complete_trio_samples = [trio.s for trio in pedigree.complete_trios()] + \
                        [trio.pat_id for trio in pedigree.complete_trios()] + \
                        [trio.mat_id for trio in pedigree.complete_trios()]
if len(complete_trio_samples)==0:
    complete_trio_samples = ['']
trio_pedigree = pedigree.filter_to(complete_trio_samples)

td = hl.trio_matrix(case_filtered_mt, trio_pedigree, complete_trios=True)

td = parent_aware_t_u_annotations_v4(td)

td = td.annotate_entries(total_t_from_parents = td.t_from_dad + td.t_from_mom,
                        total_u_from_parents = td.u_from_dad + td.u_from_mom)

# Original function for site-level
tdt_table_filtered = hl.transmission_disequilibrium_test(case_filtered_mt, trio_pedigree)

td = td.annotate_rows(tdt=tdt_table_filtered[td.row_key])

td_mt_uri = f"{prefix}.tdt.mt"
td = td.checkpoint(td_mt_uri, overwrite=True)

td = hl.read_matrix_table(td_mt_uri)

# syn_td = td.filter_rows(td.worst_csq.most_severe_consequence=='synonymous_variant')

# syn_td = syn_td.annotate_cols(total_t_from_dad=hl.agg.sum(syn_td.t_from_dad),
#                     total_u_from_dad=hl.agg.sum(syn_td.u_from_dad),
#                     total_t_from_mom=hl.agg.sum(syn_td.t_from_mom),
#                     total_u_from_mom=hl.agg.sum(syn_td.u_from_mom))

# syn_td = syn_td.annotate_rows(
#     t_total = hl.agg.sum(syn_td.t_from_dad + syn_td.t_from_mom + syn_td.t_indeterminate),
#     u_total = hl.agg.sum(syn_td.u_from_dad + syn_td.u_from_mom + syn_td.u_indeterminate)
# )


## Ultra-rare inherited ##
# inh_td = td.filter_entries((td.total_t_from_parents==1) & (td.total_u_from_parents==0))
# inh_td = inh_td.filter_rows(hl.agg.count_where((hl.is_defined(inh_td.proband_entry.GT)) |
#                   (hl.is_defined(inh_td.mother_entry.GT)) |
#                   (hl.is_defined(inh_td.father_entry.GT)))>0)

# Filter out all t=0, u=0
inh_td = td.filter_entries((td.total_t_from_parents==0) & 
                           (td.total_u_from_parents==0) &
                           (td.t_indeterminate==0) &
                           (td.u_indeterminate==0), keep=False)

# Output coding only
if coding_only:
    inh_td = inh_td.filter_rows(hl.array(coding_variants).contains(
        inh_td.worst_csq.most_severe_consequence))
inh_td_uri = f"{prefix}.ultra.rare.inherited.mt"
inh_td = inh_td.checkpoint(inh_td_uri, overwrite=True)
inh_output_uri = f"{prefix}.ultra.rare.inherited.tsv.gz"
inh_td.entries().flatten().export(inh_output_uri)

## Ultra-rare in cases/controls ##
tm = hl.trio_matrix(ultra_rare_mt, pedigree)

# Cases not in trio
non_trio_cases = ped_ht.filter((ped_ht.phenotype==2) &
             (~hl.array(complete_trio_samples).contains(ped_ht.sample_id))).sample_id.collect()
if len(non_trio_cases)==0:
    non_trio_cases = ['']
non_trio_cases_mt = case_filtered_mt.filter_cols(hl.array(non_trio_cases).contains(case_filtered_mt.s))
non_trio_cases_mt = non_trio_cases_mt.filter_entries(non_trio_cases_mt.GT.is_non_ref())
# Output coding only
if coding_only:
    non_trio_cases_mt = non_trio_cases_mt.filter_rows(hl.array(coding_variants).contains(
        non_trio_cases_mt.worst_csq.most_severe_consequence))
non_trio_cases_mt_uri = f"{prefix}.ultra.rare.non.trio.cases.mt"
non_trio_cases_mt = non_trio_cases_mt.checkpoint(non_trio_cases_mt_uri, overwrite=True)
non_trio_cases_output_uri = f"{prefix}.ultra.rare.non.trio.cases.tsv.gz"
non_trio_cases_mt.entries().flatten().export(non_trio_cases_output_uri)

# Control/unaffected samples that aren't parents in complete trios
all_parent_samples = [s for s in [trio.pat_id for trio in pedigree.trios] + \
    [trio.mat_id for trio in pedigree.trios] if s not in ['paternal_id','maternal_id',None]]
complete_trio_parent_samples = list(np.intersect1d(all_parent_samples, complete_trio_samples))
control_samples = ped_ht.filter((ped_ht.phenotype==1)).sample_id.collect()
if len(complete_trio_parent_samples)==0:
    complete_trio_parent_samples = ['']
if len(control_samples)==0:
    control_samples = ['']
control_tm = tm.filter_cols((hl.array(control_samples).contains(tm.id)) &
                           (~hl.array(complete_trio_parent_samples).contains(tm.id)))
control_tm = control_tm.filter_entries(control_tm.proband_entry.GT.is_non_ref())
# Output coding only
if coding_only:
    control_tm = control_tm.filter_rows(hl.array(coding_variants).contains(
        control_tm.worst_csq.most_severe_consequence))
control_tm_uri = f"{prefix}.ultra.rare.controls.mt"
control_tm = control_tm.checkpoint(control_tm_uri, overwrite=True)
control_output_uri = f"{prefix}.ultra.rare.controls.tsv.gz"
control_tm.entries().flatten().export(control_output_uri)