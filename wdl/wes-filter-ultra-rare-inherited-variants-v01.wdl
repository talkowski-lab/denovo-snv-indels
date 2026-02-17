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

struct QcFilters {
    Int? min_dp
    Int? max_dp
    Int? min_gq
    Int? min_pl
    Int? female_min_dp
    Int? male_auto_min_dp
    Float? het_ab_threshold
    Float? het_pab_threshold
    Float? informative_read_threshold
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
        String hail_basic_filtering_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/main/scripts/wes_denovo_basic_filtering.py"
        Float call_rate_threshold=0.8

        # step2: Variant QC filters
        QcFilters? control_qc_filters
        QcFilters? trio_case_qc_filters
        QcFilters? nontrio_case_qc_filters

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

    QcFilters qc_filters_default = object {
        min_dp: 7,
        max_dp: 1000,
        min_gq: 25,
        min_pl: 25,
        female_min_dp: 10,
        male_auto_min_dp: 10,
        het_ab_threshold: 0.25,
        het_pab_threshold: 0.000000001,
        informative_read_threshold: 0.9,
        phwe_threshold: 0.000000000001
    }

    QcFilters control_qc_filters_override = select_first([control_qc_filters, qc_filters_default])
    QcFilters trio_case_qc_filters_override = select_first([trio_case_qc_filters, qc_filters_default])
    QcFilters nontrio_case_qc_filters_override = select_first([nontrio_case_qc_filters, qc_filters_default])

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
            call step2.hailBasicFilteringRemote as step2Controls {
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
                    min_dp=select_first(
                        [
                            control_qc_filters_override.min_dp,
                            qc_filters_default.min_dp
                        ]
                    ),
                    max_dp=select_first(
                        [
                            control_qc_filters_override.max_dp,
                            qc_filters_default.max_dp
                        ]
                    ),
                    min_gq=select_first(
                        [
                            control_qc_filters_override.min_gq,
                            qc_filters_default.min_gq
                        ]
                    ),
                    min_pl=select_first(
                        [
                            control_qc_filters_override.min_pl,
                            qc_filters_default.min_pl
                        ]
                    ),
                    female_min_dp=select_first(
                        [
                            control_qc_filters_override.female_min_dp,
                            qc_filters_default.female_min_dp
                        ]
                    ),
                    male_auto_min_dp=select_first(
                        [
                            control_qc_filters_override.male_auto_min_dp,
                            qc_filters_default.male_auto_min_dp
                        ]
                    ),
                    het_ab_threshold=select_first(
                        [
                            control_qc_filters_override.het_ab_threshold,
                            qc_filters_default.het_ab_threshold
                        ]
                    ),
                    het_pab_threshold=select_first(
                        [
                            control_qc_filters_override.het_pab_threshold,
                            qc_filters_default.het_pab_threshold
                        ]
                    ),
                    informative_read_threshold=select_first(
                        [
                            control_qc_filters_override.informative_read_threshold,
                            qc_filters_default.informative_read_threshold
                        ]
                    ),
                    phwe_threshold=select_first(
                        [
                            control_qc_filters_override.phwe_threshold,
                            qc_filters_default.phwe_threshold
                        ]
                    )
            }
        }

        if (length(read_lines(subsetTrioCaseControlPed.trio_cases_ped)) > 1) {
            call step2.hailBasicFilteringRemote as step2TrioCases {
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
                    min_dp=select_first(
                        [
                            trio_case_qc_filters_override.min_dp,
                            qc_filters_default.min_dp
                        ]
                    ),
                    max_dp=select_first(
                        [
                            trio_case_qc_filters_override.max_dp,
                            qc_filters_default.max_dp
                        ]
                    ),
                    min_gq=select_first(
                        [
                            trio_case_qc_filters_override.min_gq,
                            qc_filters_default.min_gq
                        ]
                    ),
                    min_pl=select_first(
                        [
                            trio_case_qc_filters_override.min_pl,
                            qc_filters_default.min_pl
                        ]
                    ),
                    female_min_dp=select_first(
                        [
                            trio_case_qc_filters_override.female_min_dp,
                            qc_filters_default.female_min_dp
                        ]
                    ),
                    male_auto_min_dp=select_first(
                        [
                            trio_case_qc_filters_override.male_auto_min_dp,
                            qc_filters_default.male_auto_min_dp
                        ]
                    ),
                    het_ab_threshold=select_first(
                        [
                            trio_case_qc_filters_override.het_ab_threshold,
                            qc_filters_default.het_ab_threshold
                        ]
                    ),
                    het_pab_threshold=select_first(
                        [
                            trio_case_qc_filters_override.het_pab_threshold,
                            qc_filters_default.het_pab_threshold
                        ]
                    ),
                    informative_read_threshold=select_first(
                        [
                            trio_case_qc_filters_override.informative_read_threshold,
                            qc_filters_default.informative_read_threshold
                        ]
                    ),
                    phwe_threshold=select_first(
                        [
                            trio_case_qc_filters_override.phwe_threshold,
                            qc_filters_default.phwe_threshold
                        ]
                    )
            }
        }

        if (length(read_lines(subsetTrioCaseControlPed.nontrio_cases_ped)) > 1) {
            call step2.hailBasicFilteringRemote as step2NontrioCases {
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
                    min_dp=select_first(
                        [
                            nontrio_case_qc_filters_override.min_dp,
                            qc_filters_default.min_dp
                        ]
                    ),
                    max_dp=select_first(
                        [
                            nontrio_case_qc_filters_override.max_dp,
                            qc_filters_default.max_dp
                        ]
                    ),
                    min_gq=select_first(
                        [
                            nontrio_case_qc_filters_override.min_gq,
                            qc_filters_default.min_gq
                        ]
                    ),
                    min_pl=select_first(
                        [
                            nontrio_case_qc_filters_override.min_pl,
                            qc_filters_default.min_pl
                        ]
                    ),
                    female_min_dp=select_first(
                        [
                            nontrio_case_qc_filters_override.female_min_dp,
                            qc_filters_default.female_min_dp
                        ]
                    ),
                    male_auto_min_dp=select_first(
                        [
                            nontrio_case_qc_filters_override.male_auto_min_dp,
                            qc_filters_default.male_auto_min_dp
                        ]
                    ),
                    het_ab_threshold=select_first(
                        [
                            nontrio_case_qc_filters_override.het_ab_threshold,
                            qc_filters_default.het_ab_threshold
                        ]
                    ),
                    het_pab_threshold=select_first(
                        [
                            nontrio_case_qc_filters_override.het_pab_threshold,
                            qc_filters_default.het_pab_threshold
                        ]
                    ),
                    informative_read_threshold=select_first(
                        [
                            nontrio_case_qc_filters_override.informative_read_threshold,
                            qc_filters_default.informative_read_threshold
                        ]
                    ),
                    phwe_threshold=select_first(
                        [
                            nontrio_case_qc_filters_override.phwe_threshold,
                            qc_filters_default.phwe_threshold
                        ]
                    )
            }
        }

        # call step2.hailBasicFilteringRemote as step2NontrioCases {
        #     input:
        #         lcr_uri=lcr_uri,
        #         annot_mt=step1.annot_mt,
        #         input_size=getStep1MTSize.mt_size,
        #         ped_sex_qc=ped_sex_qc,
        #         bucket_id=bucket_id,
        #         cohort_prefix=cohort_prefix,
        #         hail_basic_filtering_script=hail_basic_filtering_script,
        #         call_rate_threshold=call_rate_threshold,
        #         genome_build=genome_build,
        #         hail_docker=hail_docker,
        #         # Hardcoded filters
        #         min_dp=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.min_dp,
        #                 qc_filters_default.min_dp
        #             ]
        #         ),
        #         max_dp=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.max_dp,
        #                 qc_filters_default.max_dp
        #             ]
        #         ),
        #         min_gq=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.min_gq,
        #                 qc_filters_default.min_gq
        #             ]
        #         ),
        #         min_pl=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.min_pl,
        #                 qc_filters_default.min_pl
        #             ]
        #         ),
        #         female_min_dp=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.female_min_dp,
        #                 qc_filters_default.female_min_dp
        #             ]
        #         ),
        #         male_auto_min_dp=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.male_auto_min_dp,
        #                 qc_filters_default.male_auto_min_dp
        #             ]
        #         ),
        #         het_ab_threshold=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.het_ab_threshold,
        #                 qc_filters_default.het_ab_threshold
        #             ]
        #         ),
        #         het_pab_threshold=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.het_pab_threshold,
        #                 qc_filters_default.het_pab_threshold
        #             ]
        #         ),
        #         informative_read_threshold=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.informative_read_threshold,
        #                 qc_filters_default.informative_read_threshold
        #             ]
        #         ),
        #         phwe_threshold=select_first(
        #             [
        #                 nontrio_case_qc_filters_override.phwe_threshold,
        #                 qc_filters_default.phwe_threshold
        #             ]
        #         )
        # }

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
                # bucket_id='test',
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
            <(echo -e ~{sep='\n' na_vals}) ~{ped}
            > ~{cohort_prefix}.nontrio_cases.ped

        awk -F'\t' 'NR == FNR {nas[$1]; next} (FNR == 1) || (($6 == 2) && !($3 in nas) && !($4 in nas))'
            <(echo -e ~{sep='\n' na_vals}) ~{ped}
            > ~{cohort_prefix}.trio_cases.ped
    >>>

    output {
        File controls_ped = "~{cohort_prefix}.controls.ped"
        File trio_cases_ped = "~{cohort_prefix}.trio_cases.ped"
        File nontrio_cases_ped = "~{cohort_prefix}.nontrio_cases.ped"
    }
}