task plot {

    Int loglog_ylim
    String docker

    File sumstat_file
    File conf
    Int i

    String base = basename(sumstat_file)
    String dollar="$"

    command <<<

        # Get column names from conf file
        POS=$(grep '"pos":' ${conf} | sed -n ${i}p | sed 's/.*"pos"://;s/^[^"]*"//;s/".*//')
        CHR=$(grep '"chr":' ${conf} | sed -n ${i}p | sed 's/.*"chr"://;s/^[^"]*"//;s/".*//')
        PVAL=$(grep '"pval":' ${conf} | sed -n ${i}p | sed 's/.*"pval"://;s/^[^"]*"//;s/".*//')

        mv ${sumstat_file} ${base}

        qqplot.R --file ${base} --bp_col ${dollar}{POS} --chrcol ${dollar}{CHR} --pval_col ${dollar}{PVAL} --loglog_ylim ${loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: 15*ceil(size(sumstat_file, "G")) + " GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}

workflow run_qc {

    String conf
    Array[String] sumstat_files

    Int n_studies = length(sumstat_files)

    scatter (i in range(n_studies)) {
        call plot {
            input: sumstat_file=sumstat_files[i], conf=conf, i=i+1
        }
    }
}