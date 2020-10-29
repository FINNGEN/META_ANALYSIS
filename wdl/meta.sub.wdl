task run_range {

    # localize the files, referenced locally in the conf
    Array[File] sumstat_files

    String docker
    String method
    String opts

    String pheno
    String chrom
    File conf

    command <<<

        echo "`date` FinnGen - UKBB meta-analysis"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"
        echo "method: ${method}"
        echo "options: ${opts}"
        echo "chromosome: ${chrom}"
        echo "conf: ${conf}"

        /META_ANALYSIS/scripts/meta_analysis.py ${conf} ${pheno}_chr${chrom}_meta_out.tsv ${method} ${opts} --chrom ${chrom}

        echo "`date` done"
    >>>

    output {
        File out = pheno + "_chr" + chrom + "_meta_out.tsv.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    }
}

task combine_chrom_metas {

    String pheno
    String docker
    Array[File] meta_outs

    command <<<

        echo "`date` Combining metadata files"
        echo "docker: ${docker}"
        echo "pheno: ${pheno}"

        cat <(zcat ${meta_outs[0]} | head -1) \
            <(for file in ${sep=" " meta_outs}; do
                zcat $file | tail -n +2;
            done) \
        | bgzip > ${pheno}_meta_out.tsv.gz

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ${pheno}_meta_out.tsv.gz
        echo "`date` done"
    >>>

    output {
        File meta_out = pheno + "_meta_out.tsv.gz"
        File meta_out_tbi = pheno + "_meta_out.tsv.gz.tbi"
    }
    
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    }
}


workflow run_meta {

    String pheno
    String conf
    Array[String] sumstat_files
    
    scatter (chr in range(23)) {
        call run_range {
            input: pheno=pheno, conf=conf, sumstat_files = sumstat_files, chrom=chr+1
        }
    }

    call combine_chrom_metas {
        input: pheno=pheno, meta_outs=run_range.out
    }
}