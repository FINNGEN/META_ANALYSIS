import "meta.sub.wdl" as sub

workflow meta_analysis {

    File sumstats_loc
    File pheno_confs

    Array[Array[String]] sumstat_files = read_tsv(sumstats_loc)
    Array[String] pheno_conf = read_lines(pheno_confs)

    scatter (i in range(length(pheno_conf))) {

        String pheno = basename(pheno_conf[i], ".json")

        call sub.run_meta as run_meta {
            input: pheno = pheno, conf = pheno_conf[i], sumstat_files = sumstat_files[i]
        }
    }

    output {
        Array[File] metas = run_meta.meta_out
        Array[File] filtered_metas = run_meta.filtered_meta_out
        Array[Array[File]] meta_pngs = run_meta.pngs
        Array[Array[File]] meta_lambdas = run_meta.lambdas
    }
}