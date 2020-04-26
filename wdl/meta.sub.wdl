task run_range {

    String docker
    String pheno
    String method
    File conf
    # not used but needed to localize files - locally referenced to in the conf file
    Array[File] summary_stats
    String opts
    String chrom
    Int n_studies = length(summary_stats)

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - run meta"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"
        echo "method: ${method}"
        echo "options: ${opts}"
        echo "chrom: ${chrom}"
        echo "conf: ${conf}"
        echo "n studies: ${n_studies}"
        printf "${sep='\n' summary_stats}\n"

        /META_ANALYSIS/scripts/meta_analysis.py ${opts} --chrom ${chrom} ${conf} "${pheno}_${method}_chr${chrom}_meta" ${method} && \

        echo "`date` done"

    >>>

    output {
        File out = pheno + "_" + method + "_chr" + chrom + "_meta.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

task gather {

    String docker
    Array[File] meta_stats
    String pheno
    String method
    File conf
    Array[String] summary_stats
    String opts
    Float p_thresh
    Int n_studies = length(summary_stats)
    Int n_pieces = length(meta_stats)
    Int loglog_ylim

    command <<<

        echo "`date` COVID-19 HGI meta-analysis - gather results"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"
        echo "method: ${method}"
        echo "options: ${opts}"
        echo "conf: ${conf}"
        echo "n result pieces: ${n_pieces}"
        echo "n studies: ${n_studies}"
        printf "${sep='\n' summary_stats}\n"

        echo "`date` gathering result pieces into one"
        gunzip -c ${meta_stats[0]} | head -1 > ${pheno}_${method}_meta
        for file in ${sep=" " meta_stats}; do
            gunzip -c $file | tail -n+2 >> ${pheno}_${method}_meta
        done

        echo "`date` filtering p-value ${p_thresh}"
        awk '
        NR==1 {for (i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["all_${method}_meta_p"] < ${p_thresh}
        ' ${pheno}_${method}_meta | \
        bgzip > ${pheno}_${method}_meta_${p_thresh}.txt.gz

        awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) a[$i]=i;}
        {print $a["#CHR"],$a["POS"],$a["all_${method}_meta_p"]}
        ' ${pheno}_${method}_meta \
        > ${pheno}_${method}_meta_p

        echo "`date` bgzipping and tabixing"
        bgzip ${pheno}_${method}_meta
        tabix -s 1 -b 2 -e 2 ${pheno}_${method}_meta.gz

        echo "`date` plotting qq and manhattan"
        qqplot.R --file ${pheno}_${method}_meta_p --bp_col "POS" --chrcol "#CHR" --pval_col "all_${method}_meta_p" --loglog_ylim ${loglog_ylim}

        echo "`date` done"

    >>>

    output {
        File out = pheno + "_" + method + "_meta.gz"
        File out_tbi = pheno + "_" + method + "_meta.gz.tbi"
        File sign = pheno + "_" + method + "_meta_" + p_thresh + ".txt.gz"
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow run_meta {

    String pheno
    String conf
    String method
    String opts
    Array[String] summary_stats

    scatter (chr in range(23)) {
        call run_range {
            input: pheno=pheno, method=method, opts=opts, conf=conf, summary_stats=summary_stats, chrom=chr+1
        }
    }

    call gather {
        input: pheno=pheno, method=method, opts=opts, conf=conf, summary_stats=summary_stats, meta_stats=run_range.out
    }
}
