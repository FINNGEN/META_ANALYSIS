task filter_af_info {

    String docker
    File file
    String af_col
    Float min_af
    String info_col
    Float min_info
    String outfile = basename(file, ".gz") + ".AF." + min_af + ".INFO." + min_info + ".gz"
    String dollar = "$"

    command <<<

        echo "COVID-19 HGI meta-analysis - filter sumstats" >> ${outfile}.log
        echo "${file}" >> ${outfile}.log
        echo "${af_col} ${min_af}" >> ${outfile}.log
        echo "${info_col} ${min_info}" >> ${outfile}.log
        echo "" >> ${outfile}.log

        echo "`date` original number of variants" >> ${outfile}.log
        gunzip -c ${file} | tail -n+2 | wc -l >> ${outfile}.log

        gunzip -c ${file} | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["${af_col}"]>${min_af} && (1-$a["${af_col}"])>${min_af} && $a["${info_col}"]>${min_info}' | \
        bgzip > ${outfile}

        echo "`date` new number of variants" >> ${outfile}.log
        gunzip -c ${outfile} | tail -n+2 | wc -l >> ${outfile}.log

        chr_col=${dollar}(gunzip -c ${outfile} | head -1 | tr '\t' '\n' | grep -nx "#CHR" | head -1 | cut -d ':' -f1)
        pos_col=${dollar}(gunzip -c ${outfile} | head -1 | tr '\t' '\n' | grep -nx "POS" | head -1 | cut -d ':' -f1)
        tabix -S 1 -s ${dollar}chr_col -b ${dollar}pos_col -e ${dollar}pos_col ${outfile}

        echo "`date` plotting qq and manhattan" >> ${outfile}.log
        qqplot.R --file ${outfile} --bp_col "POS" --chrcol "#CHR" --pval_col "p.value"
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
        cpu: 1
        memory: "20 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow filter_sumstats {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (file in sumstat_files) {
        call filter_af_info {
            input: file=file
        }
    }
}
