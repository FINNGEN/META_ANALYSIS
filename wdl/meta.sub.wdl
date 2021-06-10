workflow run_meta {

    String pheno
    String conf
    Array[String] sumstat_files
    
    scatter (chr in range(23)) {
        call run_range {
            input:
                pheno=pheno,
                conf=conf,
                sumstat_files = sumstat_files,
                chrom=chr+1
        }
    }

    call combine_chrom_metas {
        input:
            pheno=pheno,
            meta_outs=run_range.out
    }

    call add_rsids {
        input:
            meta_file=combine_chrom_metas.meta_out
    }

    call meta_qq {
        input:
            meta_file=combine_chrom_metas.meta_out
    }

    call post_filter {
        input:
            meta_file=add_rsids.meta_out
    }

    output {
        File meta_out = combine_chrom_metas.meta_out
        File meta_with_rsids_out = add_rsids.meta_out
        File filtered_meta_out = post_filter.filtered_meta_out
        Array[File] pngs = meta_qq.pngs
        Array[File] lambdas = meta_qq.lambdas
    }
}

# Run meta-analysis for each chromosome separately
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

        echo "`date` GWAS meta-analysis"
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
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Combine separately run meta-analysis result files
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
        memory: "2 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Add rsids
task add_rsids {

    File ref_file
    String docker

    File meta_file

    String base = basename(meta_file, ".tsv.gz")

    command <<<

        set -euxo pipefail

        echo "`date` Adding rsids"

        python3 <<EOF | bgzip > ${base}.tsv.gz

        import gzip

        fp_ref = gzip.open('${ref_file}', 'rt')
        ref_has_lines = True
        ref_chr = 1
        ref_pos = 0
        ref_line = fp_ref.readline()
        while ref_line.startswith("##"):
            ref_line = fp_ref.readline()
        if ref_line.startswith('#'):
            assert ref_line.rstrip('\r\n').split('\t') == '#CHROM POS ID REF ALT QUAL FILTER INFO'.split(), repr(ref_line)
        ref_h_idx = {h:i for i,h in enumerate(ref_line.rstrip('\r\n').split('\t'))}

        with gzip.open('${meta_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            print(header + '\trsid')
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = int(s[h_idx['#CHR']])
                pos = int(s[h_idx['POS']])
                ref = s[h_idx['REF']]
                alt = s[h_idx['ALT']]
                ref_vars = []
                while ref_has_lines and int(ref_chr) < chr or (int(ref_chr) == chr and ref_pos < pos):
                    ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chr and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False

                rsid = 'NA'
                for r in ref_vars:
                    if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                        rsid = r[ref_h_idx['ID']]
                        break

                print(line + '\t' + rsid)

        EOF

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ${base}.tsv.gz
        echo "`date` done"

    >>>

    output {
        File meta_out = base + ".tsv.gz"
        File out_tbi = base + ".tsv.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 2*ceil(size(meta_file, "G") + size(ref_file, "G")) + " SSD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }

}

# Generate qq and manhattan plots from meta-analysis results
task meta_qq {

    Int loglog_ylim
    String docker
    String pvals_to_plot

    File meta_file

    String base = basename(meta_file)

    command <<<

        set -euxo pipefail

        mv ${meta_file} ${base}

        /META_ANALYSIS/scripts/qqplot.R --file ${base} --bp_col "POS" --chrcol "#CHR" --pval_col ${pvals_to_plot} --loglog_ylim ${loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
        Array[File] lambdas = glob("*.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk " + 10*ceil(size(meta_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Filter out variants not in the left-most study (usually Finngen)
task post_filter {

    String docker

    File meta_file

    String base = basename(meta_file, ".tsv.gz")

    command <<<

        set -exo pipefail

        # Use the first '_beta' suffix column as the beta of the left-most variant. If NA --> remove variant
        zcat ${meta_file} | awk -v OFS='\t' '
        NR==1 {for(i=1;i<=NF;i++) if($i~"_beta$") {beta_col=i; break} print $0}
        (NR>1 && $beta_col != "NA") {print}
        ' | bgzip > ${base}_filtered.tsv.gz
        tabix -s 1 -b 2 -e 2 ${base}_filtered.tsv.gz

    >>>

    output {
        File filtered_meta_out = base + "_filtered.tsv.gz"
        File filtered_meta_out_tbi = base + "_filtered.tsv.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 3*ceil(size(meta_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}
