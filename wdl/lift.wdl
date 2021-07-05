workflow munge {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        call lift {
            input: sumstat_file=sumstat_file
        }
    }

    output {
        Array[File] lifted_sumstats = lift.out
    }
}

# liftover to 38 if needed
task lift {

    File sumstat_file

    String docker
    File b37_ref
    File b38_ref

    String chr_col
    String pos_col
    String ref_col
    String alt_col

    File tbi_file = sumstat_file + ".tbi"
    String base = basename(sumstat_file)

    command <<<

        set -euxo pipefail

        echo "GWAS meta-analysis - lift over sumstats if needed"
        echo "${sumstat_file}"
        echo "${b37_ref}"
        echo "${b38_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${tbi_file} ${base}.tbi

        tabix -R ${b37_ref} ${base} | wc -l > b37.txt
        tabix -R ${b38_ref} ${base} | wc -l > b38.txt

        echo "`date` `cat b37.txt` chr 21 variants build 37"
        echo "`date` `cat b38.txt` chr 21 variants build 38"

        if ((`cat b37.txt` == 0 && `cat b38.txt` == 0)); then
            echo "`date` no chr 21 variants found in either build, quitting"
            exit 1
        fi

        if ((`cat b37.txt` > `cat b38.txt`)); then
            echo "`date` lifting to build 38"
            time lift.py --info ${chr_col} ${pos_col} ${ref_col} ${alt_col} --chainfile /hg19ToHg38.over.chain.gz --temp_dir /cromwell_root/ --no_duplicates ${base}
            gunzip -c ${base}.lifted.gz | \
            cut -f2- | awk '
            BEGIN { FS=OFS="\t" }
            NR==1 {
                for (i=1;i<=NF;i++) {
                    sub("^anew_", "b37_", $i)
                    a[$i]=i;
                }
                print $0
            } NR>1 {
                temp=$a["${chr_col}"]; $a["${chr_col}"]=$a["b37_chr"]; $a["b37_chr"]=temp;
                temp=$a["${pos_col}"]; $a["${pos_col}"]=$a["b37_pos"]; $a["b37_pos"]=temp;
                temp=$a["${ref_col}"]; $a["${ref_col}"]=$a["b37_REF"]; $a["b37_REF"]=temp;
                temp=$a["${alt_col}"]; $a["${alt_col}"]=$a["b37_ALT"]; $a["b37_ALT"]=temp;
                sub("^0", "", $a["${chr_col}"]); sub("^chr", "", $a["${chr_col}"]); sub("^X", "23", $a["${chr_col}"]); sub("^Y", "24", $a["${chr_col}"]);
                if ($a["${chr_col}"] ~ /^[0-9]+$/) {
                    print $0
                }
            }' | bgzip > ${base}.GRCh38.tsv.gz

            tabix -S 1 -s 1 -b 2 -e 2 ${base}.GRCh38.tsv.gz
        else
            echo "`date` presumably already in build 38"
        fi

    >>>

    output {
        File out = base + ".GRCh38.tsv.gz"
        File out_tbi = base + ".GRCh38.tsv.gz.tbi"
        File lift_errors = "errors"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 10*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}