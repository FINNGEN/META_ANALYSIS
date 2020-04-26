import "meta.sub.wdl" as sub

workflow meta_analysis {

    File phenos
    File summary_files
    String conf_template

    Array[String] pheno = read_lines(phenos)
    Array[Array[String]] summary_stats = read_tsv(summary_files)

    scatter (i in range(length(pheno))) {
        String conf = sub(conf_template, "\\{PHENO\\}", pheno[i])
        call sub.run_meta {
            input: pheno = pheno[i], conf = conf, summary_stats = summary_stats[i]
        }
    }
}
