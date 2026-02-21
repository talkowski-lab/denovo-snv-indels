version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/lily-dev/wdl/helpers.wdl" as helpers
import "wes-denovo-step-01-remote-sharded-v01.wdl" as step1
# import "wes-denovo-step-02-remote-sharded-v01.wdl" as step2

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

struct QcFilters {
    Boolean? filter_snv_pass
    Boolean? filter_indel_pass
    Boolean? exclude_gq_filters
    Int? min_dpc
    Int? max_dpc
    Int? female_min_dp
    Int? male_auto_min_dp
    Int? qual_threshold
    Int? mq_threshold
    Float? sor_threshold_snv
    Float? readposranksum_threshold_snv
    Float? qd_threshold_snv
    Float? sor_threshold_indel
    Float? readposranksum_threshold_indel
    Float? qd_threshold_indel
    Float? max_homref_ab
    Float? min_het_snv_ab
    Float? max_het_snv_ab
    Float? min_het_indel_ab
    Float? max_het_indel_ab
    Float? het_pab_threshold
    Float? min_mean_gq
    Int? min_gq
    Int? min_pl
    Float? informative_read_threshold
    Float? call_rate_threshold
    Float? phwe_threshold
}

workflow filterUltraRareInheritedVariants {
    input {
        # step1
        Array[File] vep_vcf_files
        String mpc_ht_uri
        String gnomad_ht_uri
        String hail_annotation_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/main/scripts/wes_denovo_annotation.py"
        
        # step2
        File lcr_uri
        String hail_basic_filtering_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/lily-dev/scripts/wes_wgs_basic_filtering.py"
        Float call_rate_threshold=0.8
        
        # To help set variant QC filter defaults
        Boolean is_genome
        Boolean is_dragen

        # step2: Variant QC filters
        QcFilters control_qc_filters
        QcFilters trio_case_qc_filters
        QcFilters nontrio_case_qc_filters

        File ped_sex_qc
        String cohort_prefix
        String bucket_id
        String hail_docker
        String hail_ultra_rare_inherited_filtering_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/main/scripts/wes_ultra_rare_inherited_variants_hail.py"

        String genome_build='GRCh38'
        Float gnomad_non_neuro_af_threshold=0.001
        Float cohort_af_threshold=0.001
        Int cohort_ac_threshold=20
        Int? affected_ac_threshold
        Float? affected_af_threshold
        Boolean coding_only=true
    }

    # Define default QC filters for all cohorts
    QcFilters qc_filters_all_default = object {
        filter_snv_pass: true,
        female_min_dp: 10,
        male_auto_min_dp: 10,
        max_homref_ab: 0.05,
        het_pab_threshold: 0.000000001,
        phwe_threshold: 0.000000000001
    }
    # Define default QC filters for exome cohorts
    if (!is_genome) {
        QcFilters qc_filters_exome_default = object {
            filter_snv_pass: qc_filters_all_default.filter_snv_pass,
            female_min_dp: qc_filters_all_default.female_min_dp,
            male_auto_min_dp: qc_filters_all_default.male_auto_min_dp,
            max_homref_ab: qc_filters_all_default.max_homref_ab,
            het_pab_threshold: qc_filters_all_default.het_pab_threshold,
            phwe_threshold: qc_filters_all_default.phwe_threshold,
            filter_indel_pass: true,
            exclude_gq_filters: false,
            min_dpc: 7,
            max_dpc: 1000,
            min_het_snv_ab: 0.25,
            max_het_snv_ab: 0.75,
            min_het_indel_ab: 0.25,
            max_het_indel_ab: 0.75,
            min_gq: 25,
            min_pl: 25,
            informative_read_threshold: 0.9,
            call_rate_threshold: 0.8
        }
    }
    if (is_genome) {
        # Define default QC filters for DRAGEN genome cohorts
        # NOTE: After these relaxed filters, Eren applied hard filters after in a notebook
        if (is_dragen) {
            QcFilters qc_filters_genome_dragen_default = object {
                filter_snv_pass: qc_filters_all_default.filter_snv_pass,
                female_min_dp: qc_filters_all_default.female_min_dp,
                male_auto_min_dp: qc_filters_all_default.male_auto_min_dp,
                max_homref_ab: qc_filters_all_default.max_homref_ab,
                het_pab_threshold: qc_filters_all_default.het_pab_threshold,
                phwe_threshold: qc_filters_all_default.phwe_threshold,
                filter_indel_pass: false,
                exclude_gq_filters: true,
                min_dpc: 10,
                max_dpc: 200,
                qual_threshold: 30,
                mq_threshold: 0,
                sor_threshold_snv: 10,
                readposranksum_threshold_snv: -10,
                qd_threshold_snv: 0,
                sor_threshold_indel: 10,
                readposranksum_threshold_indel: -10,
                qd_threshold_indel: 0,
                min_het_snv_ab: 0.22,
                max_het_snv_ab: 0.78,
                min_het_indel_ab: 0.20,
                max_het_indel_ab: 0.80
            }
        }
        # Define default QC filters for non-DRAGEN genome cohorts
        if (!is_dragen) {
            QcFilters qc_filters_genome_nondragen_default = object {
                filter_snv_pass: qc_filters_all_default.filter_snv_pass,
                female_min_dp: qc_filters_all_default.female_min_dp,
                male_auto_min_dp: qc_filters_all_default.male_auto_min_dp,
                max_homref_ab: qc_filters_all_default.max_homref_ab,
                het_pab_threshold: qc_filters_all_default.het_pab_threshold,
                phwe_threshold: qc_filters_all_default.phwe_threshold,
                filter_indel_pass: false,
                exclude_gq_filters: false,
                min_dpc: 10,
                max_dpc: 200,
                qual_threshold: 150,
                mq_threshold: 50,
                sor_threshold_snv: 2.5,
                readposranksum_threshold_snv: -1.4,
                qd_threshold_snv: 3.0,
                sor_threshold_indel: 3.0,
                readposranksum_threshold_indel: -1.7,
                qd_threshold_indel: 4.0,
                min_het_snv_ab: 0.22,
                max_het_snv_ab: 0.78,
                min_het_indel_ab: 0.20,
                max_het_indel_ab: 0.80,
                min_mean_gq: 50,
                min_gq: 99
            }
        }
    }

    QcFilters qc_filters_default = select_first(
        [
            qc_filters_exome_default,
            qc_filters_genome_dragen_default,
            qc_filters_genome_nondragen_default
        ]
    )

    # Set QC filtering inputs as user-specified params if exists else
    # default params if exists else None
    # NOTE: Using `if defined() then else` here instead of select_first
    # since some of these are optional and may not be defined in
    # default or user-specified params
    QcFilters control_qc_filters_override = object {
        filter_snv_pass: (
            if defined(control_qc_filters.filter_snv_pass)
            then control_qc_filters.filter_snv_pass
            else qc_filters_default.filter_snv_pass
        ),
        filter_indel_pass: (
            if defined(control_qc_filters.filter_indel_pass)
            then control_qc_filters.filter_indel_pass
            else qc_filters_default.filter_indel_pass
        ),
        exclude_gq_filters: (
            if defined(control_qc_filters.exclude_gq_filters)
            then control_qc_filters.exclude_gq_filters
            else qc_filters_default.exclude_gq_filters
        ),
        min_dpc: (
            if defined(control_qc_filters.min_dpc)
            then control_qc_filters.min_dpc
            else qc_filters_default.min_dpc
        ),
        max_dpc: (
            if defined(control_qc_filters.max_dpc)
            then control_qc_filters.max_dpc
            else qc_filters_default.max_dpc
        ),
        female_min_dp: (
            if defined(control_qc_filters.female_min_dp)
            then control_qc_filters.female_min_dp
            else qc_filters_default.female_min_dp
        ),
        male_auto_min_dp: (
            if defined(control_qc_filters.male_auto_min_dp)
            then control_qc_filters.male_auto_min_dp
            else qc_filters_default.male_auto_min_dp
        ),
        qual_threshold: (
            if defined(control_qc_filters.qual_threshold)
            then control_qc_filters.qual_threshold
            else qc_filters_default.qual_threshold
        ),
        mq_threshold: (
            if defined(control_qc_filters.mq_threshold)
            then control_qc_filters.mq_threshold
            else qc_filters_default.mq_threshold
        ),
        sor_threshold_snv: (
            if defined(control_qc_filters.sor_threshold_snv)
            then control_qc_filters.sor_threshold_snv
            else qc_filters_default.sor_threshold_snv
        ),
        readposranksum_threshold_snv: (
            if defined(control_qc_filters.readposranksum_threshold_snv)
            then control_qc_filters.readposranksum_threshold_snv
            else qc_filters_default.readposranksum_threshold_snv
        ),
        qd_threshold_snv: (
            if defined(control_qc_filters.qd_threshold_snv)
            then control_qc_filters.qd_threshold_snv
            else qc_filters_default.qd_threshold_snv
        ),
        sor_threshold_indel: (
            if defined(control_qc_filters.sor_threshold_indel)
            then control_qc_filters.sor_threshold_indel
            else qc_filters_default.sor_threshold_indel
        ),
        readposranksum_threshold_indel: (
            if defined(control_qc_filters.readposranksum_threshold_indel)
            then control_qc_filters.readposranksum_threshold_indel
            else qc_filters_default.readposranksum_threshold_indel
        ),
        qd_threshold_indel: (
            if defined(control_qc_filters.qd_threshold_indel)
            then control_qc_filters.qd_threshold_indel
            else qc_filters_default.qd_threshold_indel
        ),
        max_homref_ab: (
            if defined(control_qc_filters.max_homref_ab)
            then control_qc_filters.max_homref_ab
            else qc_filters_default.max_homref_ab
        ),
        min_het_snv_ab: (
            if defined(control_qc_filters.min_het_snv_ab)
            then control_qc_filters.min_het_snv_ab
            else qc_filters_default.min_het_snv_ab
        ),
        max_het_snv_ab: (
            if defined(control_qc_filters.max_het_snv_ab)
            then control_qc_filters.max_het_snv_ab
            else qc_filters_default.max_het_snv_ab
        ),
        min_het_indel_ab: (
            if defined(control_qc_filters.min_het_indel_ab)
            then control_qc_filters.min_het_indel_ab
            else qc_filters_default.min_het_indel_ab
        ),
        max_het_indel_ab: (
            if defined(control_qc_filters.max_het_indel_ab)
            then control_qc_filters.max_het_indel_ab
            else qc_filters_default.max_het_indel_ab
        ),
        het_pab_threshold: (
            if defined(control_qc_filters.het_pab_threshold)
            then control_qc_filters.het_pab_threshold
            else qc_filters_default.het_pab_threshold
        ),
        min_mean_gq: (
            if defined(control_qc_filters.min_mean_gq)
            then control_qc_filters.min_mean_gq
            else qc_filters_default.min_mean_gq
        ),
        min_gq: (
            if defined(control_qc_filters.min_gq)
            then control_qc_filters.min_gq
            else qc_filters_default.min_gq
        ),
        min_pl: (
            if defined(control_qc_filters.min_pl)
            then control_qc_filters.min_pl
            else qc_filters_default.min_pl
        ),
        informative_read_threshold: (
            if defined(control_qc_filters.informative_read_threshold)
            then control_qc_filters.informative_read_threshold
            else qc_filters_default.informative_read_threshold
        ),
        call_rate_threshold: (
            if defined(control_qc_filters.call_rate_threshold)
            then control_qc_filters.call_rate_threshold
            else qc_filters_default.call_rate_threshold
        ),
        phwe_threshold: (
            if defined(control_qc_filters.phwe_threshold)
            then control_qc_filters.phwe_threshold
            else qc_filters_default.phwe_threshold
        )
    }
    QcFilters trio_case_qc_filters_override = object {
        filter_snv_pass: (
            if defined(trio_case_qc_filters.filter_snv_pass)
            then trio_case_qc_filters.filter_snv_pass
            else qc_filters_default.filter_snv_pass
        ),
        filter_indel_pass: (
            if defined(trio_case_qc_filters.filter_indel_pass)
            then trio_case_qc_filters.filter_indel_pass
            else qc_filters_default.filter_indel_pass
        ),
        exclude_gq_filters: (
            if defined(trio_case_qc_filters.exclude_gq_filters)
            then trio_case_qc_filters.exclude_gq_filters
            else qc_filters_default.exclude_gq_filters
        ),
        min_dpc: (
            if defined(trio_case_qc_filters.min_dpc)
            then trio_case_qc_filters.min_dpc
            else qc_filters_default.min_dpc
        ),
        max_dpc: (
            if defined(trio_case_qc_filters.max_dpc)
            then trio_case_qc_filters.max_dpc
            else qc_filters_default.max_dpc
        ),
        female_min_dp: (
            if defined(trio_case_qc_filters.female_min_dp)
            then trio_case_qc_filters.female_min_dp
            else qc_filters_default.female_min_dp
        ),
        male_auto_min_dp: (
            if defined(trio_case_qc_filters.male_auto_min_dp)
            then trio_case_qc_filters.male_auto_min_dp
            else qc_filters_default.male_auto_min_dp
        ),
        qual_threshold: (
            if defined(trio_case_qc_filters.qual_threshold)
            then trio_case_qc_filters.qual_threshold
            else qc_filters_default.qual_threshold
        ),
        mq_threshold: (
            if defined(trio_case_qc_filters.mq_threshold)
            then trio_case_qc_filters.mq_threshold
            else qc_filters_default.mq_threshold
        ),
        sor_threshold_snv: (
            if defined(trio_case_qc_filters.sor_threshold_snv)
            then trio_case_qc_filters.sor_threshold_snv
            else qc_filters_default.sor_threshold_snv
        ),
        readposranksum_threshold_snv: (
            if defined(trio_case_qc_filters.readposranksum_threshold_snv)
            then trio_case_qc_filters.readposranksum_threshold_snv
            else qc_filters_default.readposranksum_threshold_snv
        ),
        qd_threshold_snv: (
            if defined(trio_case_qc_filters.qd_threshold_snv)
            then trio_case_qc_filters.qd_threshold_snv
            else qc_filters_default.qd_threshold_snv
        ),
        sor_threshold_indel: (
            if defined(trio_case_qc_filters.sor_threshold_indel)
            then trio_case_qc_filters.sor_threshold_indel
            else qc_filters_default.sor_threshold_indel
        ),
        readposranksum_threshold_indel: (
            if defined(trio_case_qc_filters.readposranksum_threshold_indel)
            then trio_case_qc_filters.readposranksum_threshold_indel
            else qc_filters_default.readposranksum_threshold_indel
        ),
        qd_threshold_indel: (
            if defined(trio_case_qc_filters.qd_threshold_indel)
            then trio_case_qc_filters.qd_threshold_indel
            else qc_filters_default.qd_threshold_indel
        ),
        max_homref_ab: (
            if defined(trio_case_qc_filters.max_homref_ab)
            then trio_case_qc_filters.max_homref_ab
            else qc_filters_default.max_homref_ab
        ),
        min_het_snv_ab: (
            if defined(trio_case_qc_filters.min_het_snv_ab)
            then trio_case_qc_filters.min_het_snv_ab
            else qc_filters_default.min_het_snv_ab
        ),
        max_het_snv_ab: (
            if defined(trio_case_qc_filters.max_het_snv_ab)
            then trio_case_qc_filters.max_het_snv_ab
            else qc_filters_default.max_het_snv_ab
        ),
        min_het_indel_ab: (
            if defined(trio_case_qc_filters.min_het_indel_ab)
            then trio_case_qc_filters.min_het_indel_ab
            else qc_filters_default.min_het_indel_ab
        ),
        max_het_indel_ab: (
            if defined(trio_case_qc_filters.max_het_indel_ab)
            then trio_case_qc_filters.max_het_indel_ab
            else qc_filters_default.max_het_indel_ab
        ),
        het_pab_threshold: (
            if defined(trio_case_qc_filters.het_pab_threshold)
            then trio_case_qc_filters.het_pab_threshold
            else qc_filters_default.het_pab_threshold
        ),
        min_mean_gq: (
            if defined(trio_case_qc_filters.min_mean_gq)
            then trio_case_qc_filters.min_mean_gq
            else qc_filters_default.min_mean_gq
        ),
        min_gq: (
            if defined(trio_case_qc_filters.min_gq)
            then trio_case_qc_filters.min_gq
            else qc_filters_default.min_gq
        ),
        min_pl: (
            if defined(trio_case_qc_filters.min_pl)
            then trio_case_qc_filters.min_pl
            else qc_filters_default.min_pl
        ),
        informative_read_threshold: (
            if defined(trio_case_qc_filters.informative_read_threshold)
            then trio_case_qc_filters.informative_read_threshold
            else qc_filters_default.informative_read_threshold
        ),
        call_rate_threshold: (
            if defined(trio_case_qc_filters.call_rate_threshold)
            then trio_case_qc_filters.call_rate_threshold
            else qc_filters_default.call_rate_threshold
        ),
        phwe_threshold: (
            if defined(trio_case_qc_filters.phwe_threshold)
            then trio_case_qc_filters.phwe_threshold
            else qc_filters_default.phwe_threshold
        )
    }
    QcFilters nontrio_case_qc_filters_override = object {
        filter_snv_pass: (
            if defined(nontrio_case_qc_filters.filter_snv_pass)
            then nontrio_case_qc_filters.filter_snv_pass
            else qc_filters_default.filter_snv_pass
        ),
        filter_indel_pass: (
            if defined(nontrio_case_qc_filters.filter_indel_pass)
            then nontrio_case_qc_filters.filter_indel_pass
            else qc_filters_default.filter_indel_pass
        ),
        exclude_gq_filters: (
            if defined(nontrio_case_qc_filters.exclude_gq_filters)
            then nontrio_case_qc_filters.exclude_gq_filters
            else qc_filters_default.exclude_gq_filters
        ),
        min_dpc: (
            if defined(nontrio_case_qc_filters.min_dpc)
            then nontrio_case_qc_filters.min_dpc
            else qc_filters_default.min_dpc
        ),
        max_dpc: (
            if defined(nontrio_case_qc_filters.max_dpc)
            then nontrio_case_qc_filters.max_dpc
            else qc_filters_default.max_dpc
        ),
        female_min_dp: (
            if defined(nontrio_case_qc_filters.female_min_dp)
            then nontrio_case_qc_filters.female_min_dp
            else qc_filters_default.female_min_dp
        ),
        male_auto_min_dp: (
            if defined(nontrio_case_qc_filters.male_auto_min_dp)
            then nontrio_case_qc_filters.male_auto_min_dp
            else qc_filters_default.male_auto_min_dp
        ),
        qual_threshold: (
            if defined(nontrio_case_qc_filters.qual_threshold)
            then nontrio_case_qc_filters.qual_threshold
            else qc_filters_default.qual_threshold
        ),
        mq_threshold: (
            if defined(nontrio_case_qc_filters.mq_threshold)
            then nontrio_case_qc_filters.mq_threshold
            else qc_filters_default.mq_threshold
        ),
        sor_threshold_snv: (
            if defined(nontrio_case_qc_filters.sor_threshold_snv)
            then nontrio_case_qc_filters.sor_threshold_snv
            else qc_filters_default.sor_threshold_snv
        ),
        readposranksum_threshold_snv: (
            if defined(nontrio_case_qc_filters.readposranksum_threshold_snv)
            then nontrio_case_qc_filters.readposranksum_threshold_snv
            else qc_filters_default.readposranksum_threshold_snv
        ),
        qd_threshold_snv: (
            if defined(nontrio_case_qc_filters.qd_threshold_snv)
            then nontrio_case_qc_filters.qd_threshold_snv
            else qc_filters_default.qd_threshold_snv
        ),
        sor_threshold_indel: (
            if defined(nontrio_case_qc_filters.sor_threshold_indel)
            then nontrio_case_qc_filters.sor_threshold_indel
            else qc_filters_default.sor_threshold_indel
        ),
        readposranksum_threshold_indel: (
            if defined(nontrio_case_qc_filters.readposranksum_threshold_indel)
            then nontrio_case_qc_filters.readposranksum_threshold_indel
            else qc_filters_default.readposranksum_threshold_indel
        ),
        qd_threshold_indel: (
            if defined(nontrio_case_qc_filters.qd_threshold_indel)
            then nontrio_case_qc_filters.qd_threshold_indel
            else qc_filters_default.qd_threshold_indel
        ),
        max_homref_ab: (
            if defined(nontrio_case_qc_filters.max_homref_ab)
            then nontrio_case_qc_filters.max_homref_ab
            else qc_filters_default.max_homref_ab
        ),
        min_het_snv_ab: (
            if defined(nontrio_case_qc_filters.min_het_snv_ab)
            then nontrio_case_qc_filters.min_het_snv_ab
            else qc_filters_default.min_het_snv_ab
        ),
        max_het_snv_ab: (
            if defined(nontrio_case_qc_filters.max_het_snv_ab)
            then nontrio_case_qc_filters.max_het_snv_ab
            else qc_filters_default.max_het_snv_ab
        ),
        min_het_indel_ab: (
            if defined(nontrio_case_qc_filters.min_het_indel_ab)
            then nontrio_case_qc_filters.min_het_indel_ab
            else qc_filters_default.min_het_indel_ab
        ),
        max_het_indel_ab: (
            if defined(nontrio_case_qc_filters.max_het_indel_ab)
            then nontrio_case_qc_filters.max_het_indel_ab
            else qc_filters_default.max_het_indel_ab
        ),
        het_pab_threshold: (
            if defined(nontrio_case_qc_filters.het_pab_threshold)
            then nontrio_case_qc_filters.het_pab_threshold
            else qc_filters_default.het_pab_threshold
        ),
        min_mean_gq: (
            if defined(nontrio_case_qc_filters.min_mean_gq)
            then nontrio_case_qc_filters.min_mean_gq
            else qc_filters_default.min_mean_gq
        ),
        min_gq: (
            if defined(nontrio_case_qc_filters.min_gq)
            then nontrio_case_qc_filters.min_gq
            else qc_filters_default.min_gq
        ),
        min_pl: (
            if defined(nontrio_case_qc_filters.min_pl)
            then nontrio_case_qc_filters.min_pl
            else qc_filters_default.min_pl
        ),
        informative_read_threshold: (
            if defined(nontrio_case_qc_filters.informative_read_threshold)
            then nontrio_case_qc_filters.informative_read_threshold
            else qc_filters_default.informative_read_threshold
        ),
        call_rate_threshold: (
            if defined(nontrio_case_qc_filters.call_rate_threshold)
            then nontrio_case_qc_filters.call_rate_threshold
            else qc_filters_default.call_rate_threshold
        ),
        phwe_threshold: (
            if defined(nontrio_case_qc_filters.phwe_threshold)
            then nontrio_case_qc_filters.phwe_threshold
            else qc_filters_default.phwe_threshold
        )
    }

    # Split PED into parents, trio cases, non-trio cases to apply separate filters
    call subsetTrioCaseControlPed {
        input:
            ped=ped_sex_qc,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    scatter (vcf_file in vep_vcf_files) {
        call step1.hailAnnotateRemote as step1 {
            input:
                mt_uri=vcf_file,
                input_size=size(vcf_file, 'GB'),
                ped_sex_qc=ped_sex_qc,
                mpc_ht_uri=mpc_ht_uri,
                gnomad_ht_uri=gnomad_ht_uri,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_annotation_script=hail_annotation_script,
                genome_build=genome_build,
                hail_docker=hail_docker
        }

        call helpers.getHailMTSize as getStep1MTSize {
            input:
                mt_uri=step1.annot_mt,
                hail_docker=hail_docker
        }

        # Split MTs by controls, trio cases, non-trio cases to apply separate filters
        if (length(read_lines(subsetTrioCaseControlPed.controls_ped)) > 1) {
            call step2HailBasicFilteringRemote as step2Controls {
                input:
                    lcr_uri=lcr_uri,
                    annot_mt=step1.annot_mt,
                    input_size=getStep1MTSize.mt_size,
                    ped_sex_qc=ped_sex_qc,
                    bucket_id=bucket_id,
                    hail_basic_filtering_script=hail_basic_filtering_script,
                    genome_build=genome_build,
                    hail_docker=hail_docker,
                    filter_snv_pass=control_qc_filters_override.filter_snv_pass,
                    filter_indel_pass=control_qc_filters_override.filter_indel_pass,
                    exclude_gq_filters=control_qc_filters_override.exclude_gq_filters,
                    min_dpc=control_qc_filters_override.min_dpc,
                    max_dpc=control_qc_filters_override.max_dpc,
                    female_min_dp=control_qc_filters_override.female_min_dp,
                    male_auto_min_dp=control_qc_filters_override.male_auto_min_dp,
                    qual_threshold=control_qc_filters_override.qual_threshold,
                    mq_threshold=control_qc_filters_override.mq_threshold,
                    sor_threshold_snv=control_qc_filters_override.sor_threshold_snv,
                    readposranksum_threshold_snv=control_qc_filters_override.readposranksum_threshold_snv,
                    qd_threshold_snv=control_qc_filters_override.qd_threshold_snv,
                    sor_threshold_indel=control_qc_filters_override.sor_threshold_indel,
                    readposranksum_threshold_indel=control_qc_filters_override.readposranksum_threshold_indel,
                    qd_threshold_indel=control_qc_filters_override.qd_threshold_indel,
                    max_homref_ab=control_qc_filters_override.max_homref_ab,
                    min_het_snv_ab=control_qc_filters_override.min_het_snv_ab,
                    max_het_snv_ab=control_qc_filters_override.max_het_snv_ab,
                    min_het_indel_ab=control_qc_filters_override.min_het_indel_ab,
                    max_het_indel_ab=control_qc_filters_override.max_het_indel_ab,
                    het_pab_threshold=control_qc_filters_override.het_pab_threshold,
                    min_mean_gq=control_qc_filters_override.min_mean_gq,
                    min_gq=control_qc_filters_override.min_gq,
                    min_pl=control_qc_filters_override.min_pl,
                    informative_read_threshold=control_qc_filters_override.informative_read_threshold,
                    call_rate_threshold=control_qc_filters_override.call_rate_threshold,
                    phwe_threshold=control_qc_filters_override.phwe_threshold
            }
        }

        if (length(read_lines(subsetTrioCaseControlPed.trio_cases_ped)) > 1) {
            call step2HailBasicFilteringRemote as step2TrioCases {
                input:                    
                    lcr_uri=lcr_uri,
                    annot_mt=step1.annot_mt,
                    input_size=getStep1MTSize.mt_size,
                    ped_sex_qc=ped_sex_qc,
                    bucket_id=bucket_id,
                    hail_basic_filtering_script=hail_basic_filtering_script,
                    genome_build=genome_build,
                    hail_docker=hail_docker,
                    filter_snv_pass=trio_case_qc_filters_override.filter_snv_pass,
                    filter_indel_pass=trio_case_qc_filters_override.filter_indel_pass,
                    exclude_gq_filters=trio_case_qc_filters_override.exclude_gq_filters,
                    min_dpc=trio_case_qc_filters_override.min_dpc,
                    max_dpc=trio_case_qc_filters_override.max_dpc,
                    female_min_dp=trio_case_qc_filters_override.female_min_dp,
                    male_auto_min_dp=trio_case_qc_filters_override.male_auto_min_dp,
                    qual_threshold=trio_case_qc_filters_override.qual_threshold,
                    mq_threshold=trio_case_qc_filters_override.mq_threshold,
                    sor_threshold_snv=trio_case_qc_filters_override.sor_threshold_snv,
                    readposranksum_threshold_snv=trio_case_qc_filters_override.readposranksum_threshold_snv,
                    qd_threshold_snv=trio_case_qc_filters_override.qd_threshold_snv,
                    sor_threshold_indel=trio_case_qc_filters_override.sor_threshold_indel,
                    readposranksum_threshold_indel=trio_case_qc_filters_override.readposranksum_threshold_indel,
                    qd_threshold_indel=trio_case_qc_filters_override.qd_threshold_indel,
                    max_homref_ab=trio_case_qc_filters_override.max_homref_ab,
                    min_het_snv_ab=trio_case_qc_filters_override.min_het_snv_ab,
                    max_het_snv_ab=trio_case_qc_filters_override.max_het_snv_ab,
                    min_het_indel_ab=trio_case_qc_filters_override.min_het_indel_ab,
                    max_het_indel_ab=trio_case_qc_filters_override.max_het_indel_ab,
                    het_pab_threshold=trio_case_qc_filters_override.het_pab_threshold,
                    min_mean_gq=trio_case_qc_filters_override.min_mean_gq,
                    min_gq=trio_case_qc_filters_override.min_gq,
                    min_pl=trio_case_qc_filters_override.min_pl,
                    informative_read_threshold=trio_case_qc_filters_override.informative_read_threshold,
                    call_rate_threshold=trio_case_qc_filters_override.call_rate_threshold,
                    phwe_threshold=trio_case_qc_filters_override.phwe_threshold
            }
        }

        if (length(read_lines(subsetTrioCaseControlPed.nontrio_cases_ped)) > 1) {
            call step2HailBasicFilteringRemote as step2NontrioCases {
                input:
                    
                    lcr_uri=lcr_uri,
                    annot_mt=step1.annot_mt,
                    input_size=getStep1MTSize.mt_size,
                    ped_sex_qc=ped_sex_qc,
                    bucket_id=bucket_id,
                    hail_basic_filtering_script=hail_basic_filtering_script,
                    genome_build=genome_build,
                    hail_docker=hail_docker,
                    filter_snv_pass=nontrio_case_qc_filters_override.filter_snv_pass,
                    filter_indel_pass=nontrio_case_qc_filters_override.filter_indel_pass,
                    exclude_gq_filters=nontrio_case_qc_filters_override.exclude_gq_filters,
                    min_dpc=nontrio_case_qc_filters_override.min_dpc,
                    max_dpc=nontrio_case_qc_filters_override.max_dpc,
                    female_min_dp=nontrio_case_qc_filters_override.female_min_dp,
                    male_auto_min_dp=nontrio_case_qc_filters_override.male_auto_min_dp,
                    qual_threshold=nontrio_case_qc_filters_override.qual_threshold,
                    mq_threshold=nontrio_case_qc_filters_override.mq_threshold,
                    sor_threshold_snv=nontrio_case_qc_filters_override.sor_threshold_snv,
                    readposranksum_threshold_snv=nontrio_case_qc_filters_override.readposranksum_threshold_snv,
                    qd_threshold_snv=nontrio_case_qc_filters_override.qd_threshold_snv,
                    sor_threshold_indel=nontrio_case_qc_filters_override.sor_threshold_indel,
                    readposranksum_threshold_indel=nontrio_case_qc_filters_override.readposranksum_threshold_indel,
                    qd_threshold_indel=nontrio_case_qc_filters_override.qd_threshold_indel,
                    max_homref_ab=nontrio_case_qc_filters_override.max_homref_ab,
                    min_het_snv_ab=nontrio_case_qc_filters_override.min_het_snv_ab,
                    max_het_snv_ab=nontrio_case_qc_filters_override.max_het_snv_ab,
                    min_het_indel_ab=nontrio_case_qc_filters_override.min_het_indel_ab,
                    max_het_indel_ab=nontrio_case_qc_filters_override.max_het_indel_ab,
                    het_pab_threshold=nontrio_case_qc_filters_override.het_pab_threshold,
                    min_mean_gq=nontrio_case_qc_filters_override.min_mean_gq,
                    min_gq=nontrio_case_qc_filters_override.min_gq,
                    min_pl=nontrio_case_qc_filters_override.min_pl,
                    informative_read_threshold=nontrio_case_qc_filters_override.informative_read_threshold,
                    call_rate_threshold=nontrio_case_qc_filters_override.call_rate_threshold,
                    phwe_threshold=nontrio_case_qc_filters_override.phwe_threshold
            }
        }

    #     # call step2.hailBasicFilteringRemote as step2NontrioCases {
    #     #     input:
    #     #         lcr_uri=lcr_uri,
    #     #         annot_mt=step1.annot_mt,
    #     #         input_size=getStep1MTSize.mt_size,
    #     #         ped_sex_qc=ped_sex_qc,
    #     #         bucket_id=bucket_id,
    #     #         cohort_prefix=cohort_prefix,
    #     #         hail_basic_filtering_script=hail_basic_filtering_script,
    #     #         call_rate_threshold=call_rate_threshold,
    #     #         genome_build=genome_build,
    #     #         hail_docker=hail_docker,
    #     #         # Hardcoded filters
    #     #         min_dp=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.min_dp,
    #     #                 qc_filters_default.min_dp
    #     #             ]
    #     #         ),
    #     #         max_dp=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.max_dp,
    #     #                 qc_filters_default.max_dp
    #     #             ]
    #     #         ),
    #     #         min_gq=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.min_gq,
    #     #                 qc_filters_default.min_gq
    #     #             ]
    #     #         ),
    #     #         min_pl=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.min_pl,
    #     #                 qc_filters_default.min_pl
    #     #             ]
    #     #         ),
    #     #         female_min_dp=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.female_min_dp,
    #     #                 qc_filters_default.female_min_dp
    #     #             ]
    #     #         ),
    #     #         male_auto_min_dp=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.male_auto_min_dp,
    #     #                 qc_filters_default.male_auto_min_dp
    #     #             ]
    #     #         ),
    #     #         het_ab_threshold=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.het_ab_threshold,
    #     #                 qc_filters_default.het_ab_threshold
    #     #             ]
    #     #         ),
    #     #         het_pab_threshold=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.het_pab_threshold,
    #     #                 qc_filters_default.het_pab_threshold
    #     #             ]
    #     #         ),
    #     #         informative_read_threshold=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.informative_read_threshold,
    #     #                 qc_filters_default.informative_read_threshold
    #     #             ]
    #     #         ),
    #     #         phwe_threshold=select_first(
    #     #             [
    #     #                 nontrio_case_qc_filters_override.phwe_threshold,
    #     #                 qc_filters_default.phwe_threshold
    #     #             ]
    #     #         )
    #     # }

        # Merge all MTs together before running inherited filtering
        call helpers.mergeMTs as mergeStep2SubsetMTs {
            input: 
                mt_uris=select_all(
                    [
                        step2Controls.filtered_mt,
                        step2TrioCases.filtered_mt,
                        step2NontrioCases.filtered_mt
                    ]
                ),
                # cohort_prefix=cohort_prefix,
                bucket_id=bucket_id,
                merged_filename=cohort_prefix + '.merged',
                join_outer=true,
                hail_docker=hail_docker
        }

        call helpers.getHailMTSize as getStep2MergedMTSize {
            input:
                mt_uri=mergeStep2SubsetMTs.merged_mt,
                hail_docker=hail_docker
        }

        # call helpers.getHailMTSize as getStep2MTSize {
        #     input:
        #         mt_uri=step2NontrioCases.filtered_mt,
        #         hail_docker=hail_docker
        # }
    
        call hailUltraRareInheritedFilteringRemote {
            input:
                # filtered_mt=step2NontrioCases.filtered_mt,
                # input_size=getStep2MTSize.mt_size,
                filtered_mt=mergeStep2SubsetMTs.merged_mt,
                input_size=getStep2MergedMTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                vep_vcf_file=vcf_file,  # Fixed reference to current scattered file
                cohort_prefix=cohort_prefix,
                hail_ultra_rare_inherited_filtering_script=hail_ultra_rare_inherited_filtering_script,
                hail_docker=hail_docker,
                gnomad_non_neuro_af_threshold=gnomad_non_neuro_af_threshold,
                cohort_ac_threshold=cohort_ac_threshold,
                cohort_af_threshold=cohort_af_threshold,
                affected_ac_threshold=affected_ac_threshold,
                affected_af_threshold=affected_af_threshold,
                coding_only=coding_only
        }
    }

    call helpers.mergeResultsPython as mergeUltraRareInherited {
        input:
            tsvs=hailUltraRareInheritedFilteringRemote.ultra_rare_inherited_tsv,
            hail_docker=hail_docker,
            input_size=size(hailUltraRareInheritedFilteringRemote.ultra_rare_inherited_tsv, 'GB'),
            merged_filename=cohort_prefix+'.ultra.rare.inherited.tsv.gz'
    }

    call helpers.mergeResultsPython as mergeUltraRareNonTrioCases {
        input:
            tsvs=hailUltraRareInheritedFilteringRemote.ultra_rare_non_trio_cases_tsv,
            hail_docker=hail_docker,
            input_size=size(hailUltraRareInheritedFilteringRemote.ultra_rare_non_trio_cases_tsv, 'GB'),
            merged_filename=cohort_prefix+'.ultra.rare.non.trio.cases.tsv.gz'
    }

    call helpers.mergeResultsPython as mergeUltraRareControls {
        input:
            tsvs=hailUltraRareInheritedFilteringRemote.ultra_rare_controls_tsv,
            hail_docker=hail_docker,
            input_size=size(hailUltraRareInheritedFilteringRemote.ultra_rare_controls_tsv, 'GB'),
            merged_filename=cohort_prefix+'.ultra.rare.controls.tsv.gz'
    }

    output {
        File ultra_rare_inherited_wes_tsv = mergeUltraRareInherited.merged_tsv
        File ultra_rare_non_trio_cases_wes_tsv = mergeUltraRareNonTrioCases.merged_tsv
        File ultra_rare_controls_wes_tsv = mergeUltraRareControls.merged_tsv
    }
}

# Extract trio cases, non-trio cases, and controls into separate PEDs
task subsetTrioCaseControlPed {
    input {
        File ped
        String cohort_prefix
        Array[String] na_vals = ["0", "-9"]
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: 2,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        awk -F'\t' '(NR == 1) || ($6 == 1)' ~{ped} > ~{cohort_prefix}.controls.ped

        awk -F'\t' 'NR == FNR {nas[$1]; next} (FNR == 1) || (($6 == 2) && (($3 in nas) || ($4 in nas)))'
            <(echo -e "~{sep='\\n' na_vals}") ~{ped}
            > ~{cohort_prefix}.nontrio_cases.ped

        awk -F'\t' 'NR == FNR {nas[$1]; next} (FNR == 1) || (($6 == 2) && !($3 in nas) && !($4 in nas))'
            <(echo -e "~{sep='\\n' na_vals}") ~{ped}
            > ~{cohort_prefix}.trio_cases.ped
    >>>

    output {
        File controls_ped = "~{cohort_prefix}.controls.ped"
        File trio_cases_ped = "~{cohort_prefix}.trio_cases.ped"
        File nontrio_cases_ped = "~{cohort_prefix}.nontrio_cases.ped"
    }
}

task step2HailBasicFilteringRemote {
    input {
        File ped_sex_qc
        File lcr_uri
        Float input_size
        String annot_mt
        String bucket_id
        String hail_basic_filtering_script
        String hail_docker
        String genome_build
        # Variant filters
        Boolean? filter_snv_pass
        Boolean? filter_indel_pass
        Boolean? exclude_gq_filters
        Int? min_dpc
        Int? max_dpc
        Int? female_min_dp
        Int? male_auto_min_dp
        Int? qual_threshold
        Int? mq_threshold
        Float? sor_threshold_snv
        Float? readposranksum_threshold_snv
        Float? qd_threshold_snv
        Float? sor_threshold_indel
        Float? readposranksum_threshold_indel
        Float? qd_threshold_indel
        Float? max_homref_ab
        Float? min_het_snv_ab
        Float? max_het_snv_ab
        Float? min_het_indel_ab
        Float? max_het_indel_ab
        Float? het_pab_threshold
        Float? min_mean_gq
        Int? min_gq
        Int? min_pl
        Float? informative_read_threshold
        Float? call_rate_threshold
        Float? phwe_threshold
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -e
        curl ~{hail_basic_filtering_script} > hail_basic_filtering_script.py
        
        python3 hail_basic_filtering_script.py \
            --annot_mt ~{annot_mt} \
            --ped_uri ~{ped_sex_qc} \
            --cores ~{cpu_cores} \
            --mem ~{memory} \
            --bucket_id ~{bucket_id} \
            --lcr_uri ~{lcr_uri} \
            ~{if (defined(filter_snv_pass) && filter_snv_pass) then "--filter_snv_pass" else ""} \
            ~{
                if (defined(filter_indel_pass) && filter_indel_pass)
                then "--filter_indel_pass" else ""
            } \
            ~{
                if (defined(exclude_gq_filters) && exclude_gq_filters)
                then "--exclude_gq_filters" else ""
            } \
            ~{if defined(min_dpc) then "--min_dpc ~{min_dpc}" else ""} \
            ~{if defined(max_dpc) then "--max_dpc ~{max_dpc}" else ""} \
            ~{if defined(female_min_dp) then "--female_min_dp ~{female_min_dp}" else ""} \
            ~{if defined(male_auto_min_dp) then "--male_auto_min_dp ~{male_auto_min_dp}" else ""} \
            ~{if defined(qual_threshold) then "--qual_threshold ~{qual_threshold}" else ""} \
            ~{if defined(mq_threshold) then "--mq_threshold ~{mq_threshold}" else ""} \
            ~{if defined(sor_threshold_snv) then "--sor_threshold_snv ~{sor_threshold_snv}" else ""} \
            ~{
                if defined(readposranksum_threshold_snv)
                then "--readposranksum_threshold_snv ~{readposranksum_threshold_snv}"
                else ""
            } \
            ~{if defined(qd_threshold_snv) then "--qd_threshold_snv ~{qd_threshold_snv}" else ""} \
            ~{
                if defined(sor_threshold_indel)
                then "--sor_threshold_indel ~{sor_threshold_indel}" else ""
            } \
            ~{
                if defined(readposranksum_threshold_indel)
                then "--readposranksum_threshold_indel ~{readposranksum_threshold_indel}"
                else ""
            } \
            ~{
                if defined(qd_threshold_indel)
                then "--qd_threshold_indel ~{qd_threshold_indel}" else ""
            } \
            ~{if defined(max_homref_ab) then "--max_homref_ab ~{max_homref_ab}" else ""} \
            ~{if defined(min_het_snv_ab) then "--min_het_snv_ab ~{min_het_snv_ab}" else ""} \
            ~{if defined(max_het_snv_ab) then "--max_het_snv_ab ~{max_het_snv_ab}" else ""} \
            ~{if defined(min_het_indel_ab) then "--min_het_indel_ab ~{min_het_indel_ab}" else ""} \
            ~{if defined(max_het_indel_ab) then "--max_het_indel_ab ~{max_het_indel_ab}" else ""} \
            ~{if defined(het_pab_threshold) then "--het_pab_threshold ~{het_pab_threshold}" else ""} \
            ~{if defined(min_mean_gq) then "--min_mean_gq ~{min_mean_gq}" else ""} \
            ~{if defined(min_gq) then "--min_gq ~{min_gq}" else ""} \
            ~{if defined(min_pl) then "--min_pl ~{min_pl}" else ""} \
            ~{
                if defined(informative_read_threshold)
                then "--informative_read_threshold ~{informative_read_threshold}" else ""
            } \
            ~{
                if defined(call_rate_threshold)
                then "--call_rate_threshold ~{call_rate_threshold}" else ""
            } \
            ~{if defined(phwe_threshold) then "--phwe_threshold ~{phwe_threshold}" else ""}
    >>>

    String prefix = basename(annot_mt, "_wes_denovo_annot.mt")
    output {
        String filtered_mt = read_lines("mt_uri.txt")[0]
        File post_filter_sample_qc_info = "~{prefix}_final_annot_post_filter_qc_info.txt"
    }
}

# (hailUltraRareInheritedFilteringRemote task remains as provided)
task hailUltraRareInheritedFilteringRemote {
    input {
        File ped_sex_qc
        File vep_vcf_file
        Float input_size
        String filtered_mt
        String cohort_prefix
        String hail_ultra_rare_inherited_filtering_script
        String hail_docker
        Float gnomad_non_neuro_af_threshold
        Int cohort_ac_threshold
        Float cohort_af_threshold
        Int? affected_ac_threshold
        Float? affected_af_threshold
        Boolean coding_only
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        curl ~{hail_ultra_rare_inherited_filtering_script} > hail_ultra_rare_inherited_filtering_script.py
        python3 hail_ultra_rare_inherited_filtering_script.py \
            --ped-uri ~{ped_sex_qc} \
            --filt-mt-uri ~{filtered_mt} \
            --vep-vcf-uri ~{vep_vcf_file} \
            --gnomad-non-neuro-af-threshold ~{gnomad_non_neuro_af_threshold} \
            --cohort-ac-threshold ~{cohort_ac_threshold} \
            --cohort-af-threshold ~{cohort_af_threshold} \
            ~{if defined(affected_ac_threshold) then "--affected-ac-threshold ~{affected_ac_threshold}" else ""} \
            ~{if defined(affected_af_threshold) then "--affected-af-threshold ~{affected_af_threshold}" else ""} \
            ~{true='--coding-only' false='' coding_only} \
            --mem ~{memory}
    >>>

    String prefix = basename(filtered_mt, ".mt")
    output {
        File ultra_rare_inherited_tsv = "~{prefix}.ultra.rare.inherited.tsv.gz"
        File ultra_rare_non_trio_cases_tsv = "~{prefix}.ultra.rare.non.trio.cases.tsv.gz"
        File ultra_rare_controls_tsv = "~{prefix}.ultra.rare.controls.tsv.gz"
    }
}