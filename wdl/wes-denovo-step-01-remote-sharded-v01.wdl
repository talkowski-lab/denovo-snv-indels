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

workflow step1 {
    input {
        Array[String] mt_uris
        String gnomad_ht_uri
        String hail_annotation_script
        String hail_docker
        String bucket_id
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    scatter (mt_uri in mt_uris) {
        call helpers.getHailMTSize as getInputMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }
        
        # Find file extension/type of mt_uri (VCF or MT)
        String base_name = basename(mt_uri)
        String file_ext = if sub(base_name, "\\.mt$", "") != base_name then ".mt"
                        else if sub(base_name, "\\.vcf\\.gz$", "") != base_name then ".vcf.gz"
                        else if sub(base_name, "\\.vcf\\.bgz$", "") != base_name then ".vcf.bgz"
                        else ".unknown"

        call hailAnnotateRemote {
            input:
                mt_uri=mt_uri,
                input_size=getInputMTSize.mt_size,
                gnomad_ht_uri=gnomad_ht_uri,
                prefix=basename(mt_uri, file_ext),
                bucket_id=bucket_id,
                hail_annotation_script=hail_annotation_script,
                hail_docker=hail_docker,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_override
        }
    }

    output {
        Array[String] annot_mt = hailAnnotateRemote.annot_mt
        Array[File] sample_qc_info = hailAnnotateRemote.sample_qc_info
    }
}

task hailAnnotateRemote {
    input {
        Float input_size
        String mt_uri
        String bucket_id
        String prefix
        String gnomad_ht_uri
        String hail_annotation_script
        String hail_docker
        String genome_build
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
        curl ~{hail_annotation_script} > hail_annotation_script.py
        python3 hail_annotation_script.py \
            --mt-uri ~{mt_uri} \
            --prefix ~{prefix} \
            --gnomad-ht-uri ~{gnomad_ht_uri} \
            --cores ~{cpu_cores} \
            --mem ~{memory} \
            --bucket-id ~{bucket_id} \
            --genome-build ~{genome_build}
    }

    output {
        String annot_mt = read_lines('mt_uri.txt')[0]
        File sample_qc_info = read_lines('qc_out.txt')[0]
    }
}
