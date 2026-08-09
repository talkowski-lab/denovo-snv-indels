version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/helpers.wdl" as helpers

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
        Array[String] filtered_mt
        File ped_sex_qc
        File loeuf_file
        String bucket_id
        String hail_docker
        Float max_parent_ab=0.05
        Float min_child_ab=0.25
        Float min_dp_ratio=0.1
        Int min_gq=25
        Float min_p=0.05

        File? hail_denovo_filtering_script_override
    }

    scatter (mt_uri in filtered_mt) {
        call helpers.getHailMTSize as getStep2MTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
    
        call hailDenovoFilteringRemote {
            input:
                filtered_mt=mt_uri,
                input_size=getStep2MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                prefix=basename(mt_uri, "_wes_denovo_basic_filtering.mt"),
                loeuf_file=loeuf_file,
                hail_denovo_filtering_script_override=hail_denovo_filtering_script_override,
                hail_docker=hail_docker,
                max_parent_ab=max_parent_ab,
                min_child_ab=min_child_ab,
                min_dp_ratio=min_dp_ratio,
                min_gq=min_gq,
                min_p=min_p
        }
    }

    output {
        # step 3 output
        Array[File] de_novo_results_sharded = hailDenovoFilteringRemote.de_novo_results
        Array[File] de_novo_vep_sharded = hailDenovoFilteringRemote.de_novo_vep
        Array[String] de_novo_ht = hailDenovoFilteringRemote.de_novo_ht
        Array[String] tdt_mt = hailDenovoFilteringRemote.tdt_mt
        Array[String] tdt_parent_aware_mt = hailDenovoFilteringRemote.tdt_parent_aware_mt
    }
}


task hailDenovoFilteringRemote {
    input {
        File ped_sex_qc
        Float input_size
        String filtered_mt
        String bucket_id
        String prefix
        String loeuf_file
        String hail_docker
        Float max_parent_ab
        Float min_child_ab
        Float min_dp_ratio
        Int min_gq
        Float min_p

        File? hail_denovo_filtering_script_override
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
        python3 ~{default="/opt/scripts/wes_denovo_denovo_filtering.py" hail_denovo_filtering_script_override} \
            --filtered-mt ~{filtered_mt} \
            --prefix ~{prefix} \
            --ped-uri ~{ped_sex_qc} \
            --loeuf-file ~{loeuf_file} \
            --cores ~{cpu_cores} \
            --mem ~{memory} \
            --bucket-id ~{bucket_id} \
            --max-parent-ab ~{max_parent_ab} \
            --min-child-ab ~{min_child_ab} \
            --min-dp-ratio ~{min_dp_ratio} \
            --min-gq ~{min_gq} \
            --min-p ~{min_p} > stdout
    }

    output {
        File de_novo_results = "~{prefix}_wes_final_denovo.txt"
        File de_novo_vep = "~{prefix}_wes_final_denovo_vep.txt"
        String de_novo_ht = read_lines('mt_uri.txt')[0]
        String tdt_mt = read_lines('mt_uri.txt')[1]
        String tdt_parent_aware_mt = read_lines('mt_uri.txt')[2]
    }
}
