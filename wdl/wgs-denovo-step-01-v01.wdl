version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFs.wdl" as mergeVCFs

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
        String python_trio_sample_script
        String python_preprocess_script
        File lcr_uri
        File ped_sex_qc
        Array[File] annot_vcf_files
        String hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Int shards_per_chunk=10
        Int qual_threshold=150  # ~30 for DRAGEN
        Float sor_threshold_indel=3.0
        Float sor_threshold_snv=2.5
        Float readposranksum_threshold_indel=-1.7
        Float readposranksum_threshold_snv=-1.4
        Float qd_threshold_indel=4.0
        Float qd_threshold_snv=3.0
        Float mq_threshold=50
        Boolean filter_pass=true
        Boolean exclude_gq_filters=false
        Boolean sort_after_merge=false
        Boolean merge_split_vcf=false
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_merge_chunks
        RuntimeAttr? runtime_attr_merge_unchunked
    }

    call makeTrioSampleFiles {
        input:
            python_trio_sample_script=python_trio_sample_script,
            ped_sex_qc=ped_sex_qc,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    if (merge_split_vcf) {
        call mergeSplitVCF.splitFile as splitVEPFiles {
            input:
                file=write_lines(annot_vcf_files),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker
        }
        scatter (chunk_file in splitVEPFiles.chunks) {        
            call mergeVCFs.mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    sort_after_merge=sort_after_merge,
                    runtime_attr_override=runtime_attr_merge_chunk
            }
            call preprocessVCF as preprocessVCFChunk {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_sex_qc=ped_sex_qc,
                    vcf_uri=mergeChunk.merged_vcf_file,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    hail_docker=hail_docker,
                    qual_threshold=qual_threshold,
                    sor_threshold_indel=sor_threshold_indel,
                    sor_threshold_snv=sor_threshold_snv,
                    readposranksum_threshold_indel=readposranksum_threshold_indel,
                    readposranksum_threshold_snv=readposranksum_threshold_snv,
                    qd_threshold_indel=qd_threshold_indel,
                    qd_threshold_snv=qd_threshold_snv,
                    mq_threshold=mq_threshold,
                    filter_pass=filter_pass,
                    exclude_gq_filters=exclude_gq_filters,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs.mergeVCFs as mergeChunks {
            input:
                vcf_files=preprocessVCFChunk.preprocessed_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix+'.preprocessed',
                sort_after_merge=sort_after_merge,
                runtime_attr_override=runtime_attr_merge_chunks
        }
    }

    if (!merge_split_vcf) {
        scatter (vcf_uri in annot_vcf_files) {
            call preprocessVCF {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_sex_qc=ped_sex_qc,
                    vcf_uri=vcf_uri,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    hail_docker=hail_docker,
                    qual_threshold=qual_threshold,
                    sor_threshold_indel=sor_threshold_indel,
                    sor_threshold_snv=sor_threshold_snv,
                    readposranksum_threshold_indel=readposranksum_threshold_indel,
                    readposranksum_threshold_snv=readposranksum_threshold_snv,
                    qd_threshold_indel=qd_threshold_indel,
                    qd_threshold_snv=qd_threshold_snv,
                    mq_threshold=mq_threshold,
                    filter_pass=filter_pass,
                    exclude_gq_filters=exclude_gq_filters,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs.mergeVCFs as mergeVCFs {
            input:
                vcf_files=preprocessVCF.preprocessed_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix+'.preprocessed',
                sort_after_merge=sort_after_merge,
                runtime_attr_override=runtime_attr_merge_unchunked
        }
    }

    output {
        File meta_uri = makeTrioSampleFiles.meta_uri
        File trio_uri = makeTrioSampleFiles.trio_uri
        File merged_preprocessed_vcf_file = select_first([mergeChunks.merged_vcf_file, mergeVCFs.merged_vcf_file])
        File merged_preprocessed_vcf_idx = select_first([mergeChunks.merged_vcf_idx, mergeVCFs.merged_vcf_idx])
    }
}

task makeTrioSampleFiles {
    input {
        String python_trio_sample_script
        File ped_sex_qc
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{python_trio_sample_script} > python_trio_sample_script.py
    python3 python_trio_sample_script.py ~{ped_sex_qc} ~{cohort_prefix} 
    >>>
    
    output {
        File meta_uri = "~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{cohort_prefix}_trio_list.txt"
    }
}

task preprocessVCF {
    input {
        String python_preprocess_script
        File lcr_uri
        File ped_sex_qc
        File vcf_uri
        File meta_uri
        File trio_uri
        String hail_docker
        Int qual_threshold
        Float sor_threshold_indel
        Float sor_threshold_snv
        Float readposranksum_threshold_indel
        Float readposranksum_threshold_snv
        Float qd_threshold_indel
        Float qd_threshold_snv
        Float mq_threshold
        Boolean filter_pass
        Boolean exclude_gq_filters
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, "GB")
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
    String filename = basename(vcf_uri)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf_uri, ".vcf.gz") else basename(vcf_uri, ".vcf.bgz")

    String preprocessed_vcf_out = '~{prefix}.preprocessed.vcf.bgz'
    command <<<
        set -eou pipefail
        curl ~{python_preprocess_script} > python_preprocess_script.py
        python3 python_preprocess_script.py ~{lcr_uri} ~{ped_sex_qc} ~{meta_uri} ~{trio_uri} ~{vcf_uri} \
        ~{filter_pass} ~{exclude_gq_filters} ~{qual_threshold} ~{sor_threshold_indel} ~{sor_threshold_snv} \
        ~{readposranksum_threshold_indel} ~{readposranksum_threshold_snv} ~{qd_threshold_indel} ~{qd_threshold_snv} \
        ~{mq_threshold} ~{cpu_cores} ~{memory}
    >>>

    output {
        File preprocessed_vcf = preprocessed_vcf_out
        File preprocessed_vcf_idx = preprocessed_vcf_out + '.tbi'
    }
}
