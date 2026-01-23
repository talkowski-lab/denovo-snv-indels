version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers
import "wes-denovo-step-01-remote-sharded-v01.wdl" as step1
import "wes-denovo-step-02-remote-sharded-v01.wdl" as step2

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
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
        String hail_basic_filtering_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/main/scripts/wes_denovo_basic_filtering.py"
        Float call_rate_threshold=0.8

        # step2: Hardcoded filters defaults
        Int min_dp = 7
        Int max_dp = 1000
        Int min_gq = 25
        Int min_pl = 25
        Int female_min_dp = 10
        Int male_auto_min_dp = 10
        Float het_ab_threshold = 0.25
        Float het_pab_threshold = 0.000000001
        Float informative_read_threshold = 0.9
        Float phwe_threshold = 0.000000000001

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

        call step2.hailBasicFilteringRemote as step2 {
            input:
                lcr_uri=lcr_uri,
                annot_mt=step1.annot_mt,
                input_size=getStep1MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                cohort_prefix=cohort_prefix,
                hail_basic_filtering_script=hail_basic_filtering_script,
                call_rate_threshold=call_rate_threshold,
                genome_build=genome_build,
                hail_docker=hail_docker,
                # Hardcoded filters
                min_dp=min_dp,
                max_dp=max_dp,
                min_gq=min_gq,
                min_pl=min_pl,
                female_min_dp=female_min_dp,
                male_auto_min_dp=male_auto_min_dp,
                het_ab_threshold=het_ab_threshold,
                het_pab_threshold=het_pab_threshold,
                informative_read_threshold=informative_read_threshold,
                phwe_threshold=phwe_threshold
        }

        call helpers.getHailMTSize as getStep2MTSize {
            input:
                mt_uri=step2.filtered_mt,
                hail_docker=hail_docker
        }
    
        call hailUltraRareInheritedFilteringRemote {
            input:
                filtered_mt=step2.filtered_mt,
                input_size=getStep2MTSize.mt_size,
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