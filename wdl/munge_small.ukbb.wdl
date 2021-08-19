workflow munge_small {

    File sumstats_loc
    Array[File] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        call filter {
            input: sumstat_file=sumstat_file
        }
    }

    output {
        Array[File] filtered_sumstats = filter.out
        Array[File] filtered_sumstats_tbi = filter.tbi
    }
}

# Filter bad quality variants
task filter {

    File sumstat_file
    String docker

    String outfile = basename(sumstat_file, ".gz")+ ".tsv.gz"

    command <<<

        set -euxo pipefail

         
        zcat ${sumstat_file} | awk '
            BEGIN{FS="\t| "; OFS="\t"}
            NR==1 {
                for (i=1;i<=NF;i++) {
                    a[$i]=i;
                    if ($i=="pos") pos=i
                }
                print $0
            } NR>1 {
                if ($a["#chr"] ~ /^[0-9]+$/ && $a["pval_EUR"] > 0 && $a["beta_EUR"] < 1e6 && $a["beta_EUR"] > -1e6) {
                    printf $1
                    for (i=2; i<=NF; i++) {
                        if (i==pos) {
                            printf "\t%d", $i
                        } else {
                            printf "\t"$i
                        }
                    }
                    printf "\n"
                }
            }' | \
        bgzip > ${outfile}
        tabix -s 1 -b 2 -e 2 ${outfile}

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 3*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}
