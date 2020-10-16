import "meta.sub.wdl" as sub
import "qc.wdl" as qc


workflow meta_analysis {

    File sumstats_loc
    File pheno_confs

    Array[Array[String]] sumstat_files = read_tsv(sumstats_loc)
    Array[String] pheno_conf = read_lines(pheno_confs)

    scatter (i in range(length(pheno_conf))) {

        String pheno = basename(pheno_conf[i], ".json")

        call qc.run_qc {
            input: conf = pheno_conf[i], sumstat_files = sumstat_files[i]
        }

        call sub.run_meta {
            input: pheno = pheno, conf = pheno_conf[i], sumstat_files = sumstat_files[i]
        }
    }
}