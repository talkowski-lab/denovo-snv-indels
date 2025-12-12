version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers

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
        Array[String] filtered_mt
        File ped_sex_qc
        String cohort_prefix
        String hail_ultra_rare_inherited_filtering_script="https://raw.githubusercontent.com/talkowski-lab/denovo-snv-indels/refs/heads/main/scripts/wes_ultra_rare_inherited_variants_hail.py"
        String hail_docker

        Float gnomad_af_threshold=0.001
        Float cohort_af_threshold=0.001
        Int cohort_ac_threshold=20
    }

    scatter (mt_uri in filtered_mt) {
        call helpers.getHailMTSize as getHailMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
    
        call hailUltraRareInheritedFilteringRemote {
            input:
                filtered_mt=mt_uri,
                input_size=getHailMTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                cohort_prefix=cohort_prefix,
                hail_ultra_rare_inherited_filtering_script=hail_ultra_rare_inherited_filtering_script,
                hail_docker=hail_docker,
                gnomad_af_threshold=gnomad_af_threshold,
                cohort_af_threshold=cohort_af_threshold,
                cohort_ac_threshold=cohort_ac_threshold
        }
    }

    output {
        # step 3 output
        Array[File] ultra_rare_inherited_wes_tsvs = hailUltraRareInheritedFilteringRemote.output_tsv
    }
}


task hailUltraRareInheritedFilteringRemote {
    input {
        File ped_sex_qc
        Float input_size
        String filtered_mt
        String cohort_prefix
        String hail_ultra_rare_inherited_filtering_script
        String hail_docker
        Float gnomad_af_threshold
        Float cohort_af_threshold
        Int cohort_ac_threshold
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

    command {
        curl ~{hail_ultra_rare_inherited_filtering_script} > hail_ultra_rare_inherited_filtering_script.py
        python3 hail_ultra_rare_inherited_filtering_script.py \
            --ped-uri ~{ped_sex_qc} \
            --filt-mt-uri ~{filtered_mt} \
            --gnomad-af-threshold ~{gnomad_af_threshold} \
            --cohort-af-threshold ~{cohort_af_threshold} \
            --cohort-ac-threshold ~{cohort_ac_threshold} \
            --mem ~{memory}
    }

    String prefix = basename(filtered_mt, ".mt")
    output {
        File output_tsv = "~{prefix}.tdt.inherited.tsv.gz"
    }
}
