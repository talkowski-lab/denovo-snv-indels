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
        Array[File] split_trio_annot_vcfs
        String trio_denovo_docker
        Float minDQ = 2
        Boolean replace_missing_pl = true
        RuntimeAttr? runtime_attr_replace_missing_pl
        RuntimeAttr? runtime_attr_trio_denovo
        RuntimeAttr? runtime_attr_combine_vcfs
    }
    scatter (vcf_file in split_trio_annot_vcfs) {
        if (replace_missing_pl) {
            call replaceMissingPL {
                input:
                vcf_file=vcf_file,
                docker=trio_denovo_docker,
                runtime_attr_override=runtime_attr_replace_missing_pl
            }
        }
        call trio_denovo {
            input:
                ped_uri_trios=ped_uri_trios,
                vcf_file=select_first([replaceMissingPL.output_vcf, vcf_file]),
                trio_denovo_docker=trio_denovo_docker,
                minDQ=minDQ,
                runtime_attr_override=runtime_attr_trio_denovo
        }
    }
    call combineOutputVCFs {
        input:
            out_vcfs=trio_denovo.out_vcf,
            trio_denovo_docker=trio_denovo_docker,
            runtime_attr_override=runtime_attr_combine_vcfs
    }

    output {
        Array[File] trio_denovo_vcf = combineOutputVCFs.trio_denovo_vcf
    }
}

task trio_denovo {
    input {
        File ped_uri_trios
        File vcf_file
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
        cat <<EOF > getSamplePedigree.py
        import os
        import pandas as pd
        import sys

        ped_uri = sys.argv[1]
        ped = pd.read_csv(ped_uri, sep='\t', dtype={i: str for i in range(4)}).iloc[:,:6]
        ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        ped.index = ped.sample_id

        sample = sys.argv[2]

        parents = ped.loc[sample].iloc[2:4].tolist()

        ped.loc[parents+[sample]].to_csv(f"{sample}.ped", sep='\t', index=False, header=None) 
        EOF

        sample=$(basename "~{vcf_file}" '.vcf')
        sample="${sample%.filled.PL}"
        sample=$(echo "$sample" | awk -F "_trio_" '{print $2}')
        sample="${sample//_HP_VAF/}"

        python3 getSamplePedigree.py ~{ped_uri_trios} $sample
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
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(out_vcfs, "GB") 
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
        cat <<EOF > vcf_replace_missing_pl.py
        #!/usr/bin/env python3
        import argparse
        import gzip

        def open_maybe_gzip(path, mode):
            if path.endswith(".gz"):
                return gzip.open(path, mode + "t")
            return open(path, mode)

        def fill_and_reorder_PL(input_vcf, output_vcf):
            """
            For every record:
              - Move PL to be the last FORMAT field (some tools/pipelines assume this).
              - Replace a missing PL value ('.') with the 3-value missing form ('.,.,.').
            """
            with open_maybe_gzip(input_vcf, "r") as vin, open_maybe_gzip(output_vcf, "w") as vout:
                for line in vin:
                    if line.startswith("#"):
                        vout.write(line)
                        continue

                    fields = line.rstrip("\n").split("\t")
                    format_keys = fields[8].split(":")

                    if "PL" in format_keys:
                        # New FORMAT order: everything else, then PL last.
                        new_order = [k for k in format_keys if k != "PL"] + ["PL"]
                        fields[8] = ":".join(new_order)

                        for s in range(9, len(fields)):
                            sample_vals = fields[s].split(":")
                            # VCF allows trailing fields to be dropped when missing; pad them back.
                            while len(sample_vals) < len(format_keys):
                                sample_vals.append(".")

                            val_by_key = dict(zip(format_keys, sample_vals))

                            pl_val = val_by_key.get("PL", ".")
                            if pl_val == ".":
                                pl_val = ".,.,."
                            val_by_key["PL"] = pl_val

                            fields[s] = ":".join(val_by_key[k] for k in new_order)

                    vout.write("\t".join(fields) + "\n")

            print(f"Done! PL moved to last FORMAT field, missing PL filled as .,.,. Output: {output_vcf}")


        if __name__ == "__main__":
            parser = argparse.ArgumentParser(
                description="Reorder FORMAT fields so PL is always last, and fill missing PL with .,.,."
            )
            parser.add_argument("input_vcf", help="Input VCF file")
            parser.add_argument("output_vcf", help="Output VCF file")

            args = parser.parse_args()

            fill_and_reorder_PL(args.input_vcf, args.output_vcf)
        EOF

        python3 vcf_replace_missing_pl.py "~{vcf_file}" "~{output_filename}"
    >>>

    output {
        File output_vcf = output_filename 
    }
}
