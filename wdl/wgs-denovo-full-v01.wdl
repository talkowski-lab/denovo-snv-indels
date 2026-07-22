version 1.0

import "wgs-denovo-step-01-v01.wdl" as step1
import "wgs-denovo-step-02-v01.wdl" as step2
import "wgs-denovo-step-03-v01.wdl" as step3
import "wgs-denovo-step-04-v01.wdl" as step4
import "wgs-denovo-step-05-v01.wdl" as step5
import "wgs-denovo-step-06-v01.wdl" as step6
import "wgs-denovo-step-07-pu-only.wdl" as step7
import "filterUltraRareInheritedVariantsHail.wdl" as filterUltraRareInheritedVariantsHail
import "filterUltraRareParentsVariantsHail.wdl" as filterUltraRareParentsVariantsHail 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow wgs_denovo_full {
    input {
        File lcr_uri
        File ped_sex_qc
        File relatedness_qc
        File hg38_reference
        File hg38_reference_fai
        File hg38_reference_dict
        File repetitive_regions_bed

        Array[File] annot_vcf_files
        String cohort_prefix
        String sv_base_mini_docker
        String trio_denovo_docker
        String hail_docker
        String jvarkit_docker
        String sample_column
        Int batch_size=10

        # Note: only scripts that are shared between steps/tasks are input at the top-level
        # All other scripts and runtime_attr inputs are added at the sub-workflow/task level
        String python_trio_sample_script
        String prioritize_csq_script

        Boolean filter_pass=true
        Boolean exclude_gq_filters=false
        Boolean merge_split_vcf=false
        String genome_build='GRCh38'
        Int shards_per_chunk=10

        # Note: only filters that are shared between steps/tasks are input at the top-level
        # The following filters should be applied uniformly across step1, ultra-rare inherited, ultra-rare parents filtering
        # Other filters are meant to be modifiable based on filtering type, so left out as top-level inputs
        Int qual_threshold=150
        Float sor_threshold_indel=3.0
        Float sor_threshold_snv=2.5
        Float readposranksum_threshold_indel=-1.7
        Float readposranksum_threshold_snv=-1.4
        Float qd_threshold_indel=4.0
        Float qd_threshold_snv=3.0
        Float mq_threshold=50
        Float minDQ=2
    }

    call step1.step1 as step1 {
        input:
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            annot_vcf_files=annot_vcf_files,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
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
            merge_split_vcf=merge_split_vcf,
            shards_per_chunk=shards_per_chunk
    }

    call step2.step2 as step2 {
        input:
            merged_preprocessed_vcf_file=step1.merged_preprocessed_vcf_file,
            relatedness_qc=relatedness_qc,
            ped_sex_qc=ped_sex_qc,
            hail_docker=hail_docker
    }

    call step3.step3 as step3 {
        input:
            ped_sex_qc=ped_sex_qc,
            merged_preprocessed_vcf_file_filtered=step2.merged_preprocessed_vcf_file_filtered,
            hail_docker=hail_docker,
            cohort_prefix=cohort_prefix,
            trio_denovo_docker=trio_denovo_docker,
            batch_size=batch_size,
            hg38_reference=hg38_reference,
            hg38_reference_fai=hg38_reference_fai,
            hg38_reference_dict=hg38_reference_dict,
            jvarkit_docker=jvarkit_docker
    }

    call step4.step4 as step4 {
        input:
            ped_uri_trios=step3.ped_uri_trios,
            split_trio_annot_vcfs=step3.split_trio_annot_vcfs,
            trio_denovo_docker=trio_denovo_docker,
            minDQ=minDQ
    }

    call step5.step5 as step5 {
        input:
            ped_sex_qc=ped_sex_qc,
            split_trio_annot_vcfs=step3.split_trio_annot_vcfs,
            trio_denovo_vcf=step4.trio_denovo_vcf,
            trio_denovo_docker=trio_denovo_docker,
            cohort_prefix=cohort_prefix
    }

    call step6.step6 as step6 {
        input:
            vcf_metrics_tsv=step5.vcf_metrics_tsv,
            annot_vcf_files=annot_vcf_files,
            hail_docker=hail_docker,
            sample_column=sample_column,
            genome_build=genome_build,
            prioritize_csq_script=prioritize_csq_script
    }

    call filterUltraRareInheritedVariantsHail.filterUltraRareInheritedVariantsHail as filterUltraRareInheritedVariantsHail {
        input:
            annot_vcf_files=annot_vcf_files,
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            python_trio_sample_script=python_trio_sample_script,
            vcf_metrics_tsv_final=step6.vcf_metrics_tsv_final,
            hg38_reference=hg38_reference,
            hg38_reference_dict=hg38_reference_dict,
            hg38_reference_fai=hg38_reference_fai,
            jvarkit_docker=jvarkit_docker,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            qual_threshold=qual_threshold,
            sor_threshold_indel=sor_threshold_indel,
            sor_threshold_snv=sor_threshold_snv,
            readposranksum_threshold_indel=readposranksum_threshold_indel,
            readposranksum_threshold_snv=readposranksum_threshold_snv,
            qd_threshold_indel=qd_threshold_indel,
            qd_threshold_snv=qd_threshold_snv,
            mq_threshold=mq_threshold,
            prioritize_gnomad=false,
            prioritize_csq_script=prioritize_csq_script
    }

    call filterUltraRareParentsVariantsHail.filterUltraRareParentsVariantsHail as filterUltraRareParentsVariantsHail {
        input:
            annot_vcf_files=annot_vcf_files,
            lcr_uri=lcr_uri,
            ped_sex_qc=ped_sex_qc,
            python_trio_sample_script=python_trio_sample_script,
            vcf_metrics_tsv_final=step6.vcf_metrics_tsv_final,
            hg38_reference=hg38_reference,
            hg38_reference_dict=hg38_reference_dict,
            hg38_reference_fai=hg38_reference_fai,
            jvarkit_docker=jvarkit_docker,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            qual_threshold=qual_threshold,
            sor_threshold_indel=sor_threshold_indel,
            sor_threshold_snv=sor_threshold_snv,
            readposranksum_threshold_indel=readposranksum_threshold_indel,
            readposranksum_threshold_snv=readposranksum_threshold_snv,
            qd_threshold_indel=qd_threshold_indel,
            qd_threshold_snv=qd_threshold_snv,
            mq_threshold=mq_threshold,
            prioritize_gnomad=true,
            prioritize_csq_script=prioritize_csq_script
    }

    call step7.step7 as step7 {
        input:
            downsampled_ultra_rare_inherited=filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_Indel,
            downsampled_ultra_rare_parents=filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_Indel,
            vcf_metrics_tsv_final=step6.vcf_metrics_tsv_final,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            repetitive_regions_bed=repetitive_regions_bed,
            var_type="Indel"  # Run PU model only on Indels by default
    }

    output {
        File meta_uri = step1.meta_uri
        File trio_uri = step1.trio_uri
        File merged_preprocessed_vcf_file = step1.merged_preprocessed_vcf_file
        File merged_preprocessed_vcf_idx = step1.merged_preprocessed_vcf_idx
        File merged_preprocessed_vcf_file_filtered = step2.merged_preprocessed_vcf_file_filtered
        File merged_preprocessed_sample_qc = step2.merged_preprocessed_sample_qc
        File ped_uri_trios = step3.ped_uri_trios
        Array[File] split_trio_vcfs = step3.split_trio_vcfs
        Array[File] split_trio_annot_vcfs = step3.split_trio_annot_vcfs
        Array[File] trio_denovo_vcf = step4.trio_denovo_vcf
        File vcf_metrics_tsv = step5.vcf_metrics_tsv
        File vcf_metrics_tsv_prior_csq = step6.vcf_metrics_tsv_prior_csq
        File vcf_metrics_tsv_final = step6.vcf_metrics_tsv_final
        File vcf_metrics_tsv_final_pu = step7.vcf_metrics_tsv_final_pu
        
        File ultra_rare_inherited_tsv = filterUltraRareInheritedVariantsHail.ultra_rare_inherited_tsv
        File downsampled_ultra_rare_inherited_SNV = filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_SNV
        File downsampled_ultra_rare_inherited_Indel = filterUltraRareInheritedVariantsHail.downsampled_ultra_rare_inherited_Indel

        File ultra_rare_parents_tsv = filterUltraRareParentsVariantsHail.ultra_rare_parents_tsv
        File downsampled_ultra_rare_parents_SNV = filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_SNV
        File downsampled_ultra_rare_parents_Indel = filterUltraRareParentsVariantsHail.downsampled_ultra_rare_parents_Indel
    }    
}
