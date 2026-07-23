version 1.0 

import "annotateHPandVAF.wdl" as annotateHPandVAF

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step3 {
    input {
        File ped_sex_qc
        File merged_preprocessed_vcf_file_filtered
        String hail_docker
        String cohort_prefix
        String trio_denovo_docker
        Int batch_size=10

        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        String jvarkit_docker

        Boolean subset_ped=true
        RuntimeAttr? runtime_attr_subset_ped
        RuntimeAttr? runtime_attr_uber_split
        RuntimeAttr? runtime_attr_annotate_hp_vaf
        RuntimeAttr? runtime_attr_combine_vcfs

    }

    String stats_file = cohort_prefix + "_stats.txt"

    if (subset_ped) {
        call subsetPed {
            input:
                ped_sex_qc=ped_sex_qc,
                vcf_file=merged_preprocessed_vcf_file_filtered,
                trio_denovo_docker=trio_denovo_docker,
                runtime_attr_override=runtime_attr_subset_ped
        }
    }

    File new_ped_sex_qc = select_first([subsetPed.new_ped_sex_qc, ped_sex_qc])

    call uberSplit_v3 {
        input:
            ped_sex_qc=new_ped_sex_qc,
            vcf_file=merged_preprocessed_vcf_file_filtered,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            stats_file=stats_file,
            batch_size=batch_size,
            runtime_attr_override=runtime_attr_uber_split
    }

    call annotateHPandVAF.annotateHPandVAF as annotateHPandVAF {
        input:
            split_trio_vcfs=uberSplit_v3.split_trio_vcfs,
            vep_vcf_file=merged_preprocessed_vcf_file_filtered,
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            jvarkit_docker=jvarkit_docker,
            runtime_attr_annotate_hp_vaf=runtime_attr_annotate_hp_vaf,
            runtime_attr_combine_vcfs=runtime_attr_combine_vcfs
    }

    output {
        File ped_uri_trios = new_ped_sex_qc
        Array[File] split_trio_vcfs = uberSplit_v3.split_trio_vcfs
        Array[File] split_trio_annot_vcfs = annotateHPandVAF.split_trio_annot_vcfs
        File stats_files = uberSplit_v3.stats_file_out
    }
}

task subsetPed {
    input {
        File ped_sex_qc
        File vcf_file
        String trio_denovo_docker
        
        File? subset_ped_script_override
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB") 
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: trio_denovo_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        bcftools query -l ~{vcf_file} > samples.txt
        python3 ~{default="/opt/scripts/subset_ped.py" subset_ped_script_override} \
            samples.txt ~{ped_sex_qc} > stdout
    >>>

    output {
        File new_ped_sex_qc = basename(ped_sex_qc, '.ped')+'_subset.ped'
    }
}

task uberSplit_v3 {
    input {
        File ped_sex_qc
        File vcf_file
        String hail_docker
        String cohort_prefix
        String stats_file
        Int batch_size
        
        File? uberSplit_v3_script_override       
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, "GB")
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command {
        set -eou pipefail
        mkdir -p ~{cohort_prefix}
        python3 ~{default="/opt/scripts/uberSplit_v3.py" uberSplit_v3_script_override} \
            ~{ped_sex_qc} ~{vcf_file} ~{cohort_prefix} ~{stats_file} ~{batch_size}
    }

    output {
        Array[File] split_trio_vcfs = glob(cohort_prefix + "/*")
        File stats_file_out = stats_file
    }
}