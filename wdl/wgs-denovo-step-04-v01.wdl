version 1.0 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step4 {
    input {
        File ped_uri_trios
        Array[File] split_trio_vcfs
        String get_sample_pedigree_script
        String trio_denovo_docker
        Float minDQ
        Boolean replace_missing_pl = false
        RuntimeAttr? runtime_attr_trio_denovo
    }
    scatter (vcf_file in split_trio_vcfs) {
        if (replace_missing_pl) {
            call replaceMissingPL {
                input:
                vcf_file=vcf_file,
                docker=trio_denovo_docker,
                runtime_attr_override=runtime_attr_trio_denovo
            }
        }
        call trio_denovo {
            input:
                ped_uri_trios=ped_uri_trios,
                vcf_file=select_first([replaceMissingPL.output_vcf, vcf_file]),
                get_sample_pedigree_script=get_sample_pedigree_script,
                trio_denovo_docker=trio_denovo_docker,
                minDQ=minDQ,
                runtime_attr_override=runtime_attr_trio_denovo
        }
    }
    call combineOutputVCFs {
        input:
            out_vcfs=trio_denovo.out_vcf,
            trio_denovo_docker=trio_denovo_docker
    }

    output {
        Array[File] trio_denovo_vcf = combineOutputVCFs.trio_denovo_vcf
    }
}

task trio_denovo {
    input {
        File ped_uri_trios
        File vcf_file
        String get_sample_pedigree_script
        String trio_denovo_docker
        Float minDQ
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
        set -eou pipefail
        sample=$(basename "~{vcf_file}" '.vcf')
        sample="${sample%.filled.PL}"
        sample=$(echo "$sample" | awk -F "_trio_" '{print $2}')
        sample="${sample//_HP_VAF/}"
        curl ~{get_sample_pedigree_script} > get_sample_pedigree_script.py
        python3 get_sample_pedigree_script.py ~{ped_uri_trios} $sample
        /src/wgs_denovo/triodenovo/triodenovo-fix/src/triodenovo --ped "$sample".ped \
            --in_vcf "~{vcf_file}" \
            --out_vcf "~{basename(vcf_file, '.vcf') + '.denovos.vcf'}" \
            --minDQ ~{minDQ}
        bgzip "~{basename(vcf_file, '.vcf') + '.denovos.vcf'}"
    >>>

    output {
        File out_vcf = basename(vcf_file, '.vcf') + '.denovos.vcf.gz'
    }
}

task combineOutputVCFs {
    input {
        Array[File] out_vcfs
        String trio_denovo_docker
    }

    runtime {
        docker: trio_denovo_docker
    }

    String concat_files_string = "\\"~{sep='" "' out_vcfs}\\""

    command <<<
        set -eou pipefail
        mkdir -p tmp_out_vcfs

        # write the list of files to a temp text file
        VCFS="~{write_lines(out_vcfs)}"
        cat $VCFS > vcfs_list.txt

        # safely move each file listed in the text file
        while IFS= read -r f; do
            mv "$f" tmp_out_vcfs/
        done < "vcfs_list.txt"
    >>>

    output {
        Array[File] trio_denovo_vcf = glob('tmp_out_vcfs/*')
    }
}

task replaceMissingPL {
    input {
        File vcf_file
        String docker
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    String output_filename = basename(vcf_file, '.vcf') + '.filled.PL.vcf'
    command <<<
        set -eou pipefail
        bcftools +setGT "~{vcf_file}" -Ov -o "~{output_filename}" -- \
            -t q \
            -n '0,0,0' \
            -i 'FMT/PL="."'
    >>>

    output {
        File output_vcf = output_filename 
    }
}
