version 1.0

import "wes-denovo-step-01-remote-sharded-v01.wdl" as step1
import "wes-denovo-step-02-remote-sharded-v01.wdl" as step2
import "wes-denovo-step-03-remote-sharded-v01.wdl" as step3
import "wes-denovo-step-04-remote-sharded-v01.wdl" as step4
import "wes-denovo-step-05-remote-sharded-v01.wdl" as step5
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/helpers.wdl" as helpers
import "prioritizeCSQ.wdl" as prioritizeCSQ_og

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow hailDenovoWES {
    input {
        Array[String] mt_uris
        File lcr_uri
        File ped_sex_qc
        File loeuf_file
        File eval_regions

        String bucket_id
        String gnomad_ht_uri
        String cohort_prefix
        String genome_build

        String hail_docker

        # step2 filters
        Float call_rate_threshold=0.8
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

        # step3 filters
        Float max_parent_ab=0.05
        Float min_child_ab=0.25
        Float min_dp_ratio=0.1
        Int min_gq=25
        Float min_p=0.05
    }

    scatter (mt_uri in mt_uris) {
        # Find file extension/type of mt_uri (VCF or MT)
        String base_name = basename(mt_uri)
        String file_ext = if sub(base_name, "\\.mt$", "") != base_name then ".mt"
                        else if sub(base_name, "\\.vcf\\.gz$", "") != base_name then ".vcf.gz"
                        else if sub(base_name, "\\.vcf\\.bgz$", "") != base_name then ".vcf.bgz"
                        else ".unknown"
        String mt_uri_prefix = basename(mt_uri, file_ext)
        
        call helpers.getHailMTSize as getInputMTSize {
            input:
                mt_uri=mt_uri,
                hail_docker=hail_docker
        }

        call step1.hailAnnotateRemote as step1 {
            input:
                mt_uri=mt_uri,
                prefix=mt_uri_prefix,
                input_size=getInputMTSize.mt_size,
                gnomad_ht_uri=gnomad_ht_uri,
                bucket_id=bucket_id,
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
                prefix=mt_uri_prefix,
                input_size=getStep1MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                genome_build=genome_build,
                hail_docker=hail_docker,
                # Passing parameters to task
                call_rate_threshold=call_rate_threshold,
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
    
        call step3.hailDenovoFilteringRemote as step3 {
            input:
                filtered_mt=step2.filtered_mt,
                prefix=mt_uri_prefix,
                input_size=getStep2MTSize.mt_size,
                ped_sex_qc=ped_sex_qc,
                bucket_id=bucket_id,
                loeuf_file=loeuf_file,
                hail_docker=hail_docker,
                max_parent_ab=max_parent_ab,
                min_child_ab=min_child_ab,
                min_dp_ratio=min_dp_ratio,
                min_gq=min_gq,
                min_p=min_p
        }
    }

    call step4.step4 as step4 {
        input:
        de_novo_results_sharded=step3.de_novo_results, 
        de_novo_vep_sharded=step3.de_novo_vep,
        vep_vcf_file=mt_uris[0],
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    call step5.step5 as step5 {
        input:
        de_novo_merged=step4.de_novo_merged,
        eval_regions=eval_regions,
        cohort_prefix=cohort_prefix,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    output {
        # step 1 output
        Array[String] annot_mt = step1.annot_mt
        Array[File] sample_qc_info = step1.sample_qc_info
        # step 2 output
        Array[String] filtered_mt = step2.filtered_mt
        Array[File] post_filter_sample_qc_info = step2.post_filter_sample_qc_info
        # step 3 output
        Array[String] de_novo_ht = step3.de_novo_ht
        Array[String] tdt_mt = step3.tdt_mt
        Array[String] tdt_parent_aware_mt = step3.tdt_parent_aware_mt
        # step4 output
        File de_novo_results = step4.de_novo_results
        File de_novo_vep = step4.de_novo_vep
        File de_novo_merged = step4.de_novo_merged
        # step5 output
        File de_novo_final = step5.de_novo_final
    }
}