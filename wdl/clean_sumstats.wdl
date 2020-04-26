task clean {

    String docker
    File sumstat_file
    String outfile = sub(basename(sumstat_file, ".gz"), "\\.bgz$", "") + ".munged.gz"
    String dollar = "$"

    command <<<

        echo "COVID-19 HGI meta-analysis - clean sumstats" >> ${outfile}.log
        echo "${sumstat_file}" >> ${outfile}.log
        echo "" >> ${outfile}.log

        echo "`date` original number of variants" >> ${outfile}.log
        gunzip -c ${sumstat_file} | tail -n+2 | wc -l >> ${outfile}.log

        chr_col=${dollar}(gunzip -c ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "CHR" | head -1 | cut -d ':' -f1)
        pos_col=${dollar}(gunzip -c ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "POS" | head -1 | cut -d ':' -f1)
        printf "`date` col CHR "${dollar}{chr_col}" col POS "${dollar}{pos_col}"\n" >> ${outfile}.log

        gunzip -c ${sumstat_file} | awk ' \
          BEGIN{FS="\t| "; OFS="\t"}
          NR==1 {
              for (i=1;i<=NF;i++) { sub("^CHR", "#CHR", $i); a[$i]=i; if ($i=="POS") pos=i }
              gsub("Pvalue", "p.value", $0);
              print $0
          } NR>1 {
              sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]);
              if ($a["#CHR"] ~ /^[0-9]+$/) {
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
        sort -k${dollar}chr_col,${dollar}{chr_col}g -k${dollar}pos_col,${dollar}{pos_col}g -u | \
        bgzip > ${outfile}

        echo "`date` new number of variants" >> ${outfile}.log
        gunzip -c ${outfile} | tail -n+2 | wc -l >> ${outfile}.log
        echo "`date` headers" >> ${outfile}.log
        gunzip -c ${outfile} | head -1 | tr '\t' '\n' >> ${outfile}.log

        tabix -S 1 -s ${dollar}chr_col -b ${dollar}pos_col -e ${dollar}pos_col ${outfile}

        gunzip -c ${outfile} | tail -n+2 | cut -f ${dollar}chr_col | uniq > chr.tmp
        echo "`date` ${dollar}(wc -l chr.tmp | cut -d' ' -f1) chromosomes" >> ${outfile}.log
        cat chr.tmp >> ${outfile}.log

        echo "`date` plotting qq and manhattan" >> ${outfile}.log
        qqplot.R --file ${outfile} --bp_col "POS" --chrcol "#CHR" --pval_col "p.value"

        echo "`date` unique number of fields" >> ${outfile}.log
        gunzip -c ${outfile} | awk 'BEGIN{FS="\t"} {print NF}' | sort -u > n.tmp
        cat n.tmp >> ${outfile}.log
        if [ $(wc -l n.tmp | cut -d' ' -f1) != 1 ]; then echo "file not square"; exit 1; fi
        if [ $(wc -l chr.tmp | cut -d' ' -f1) -lt 22 ]; then echo "less than 22 chromosomes"; exit 1; fi

        echo "`date` done" >> ${outfile}.log

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
        File log = outfile + ".log"
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: 15*ceil(size(sumstat_file, "G")) + " GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow clean_sumstats {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        call clean {
            input: sumstat_file=sumstat_file
        }
    }
}
