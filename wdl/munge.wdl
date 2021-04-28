workflow munge {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        call clean_filter {
            input: sumstat_file=sumstat_file
        }
        call lift {
            input: sumstat_file=clean_filter.out
        }
        call harmonize {
            input: sumstat_file=lift.out
        }
        call plot {
            input: sumstat_file=harmonize.out
        }
    }

    output {
        Array[File] filtered_sumstats = clean_filter.out
        Array[File] lifted_sumstats = lift.out
        Array[File] harmonized_sumstats = harmonize.out
        Array[Array[File]] plots = plot.pngs
    }
}

# Filter bad quality variants
task clean_filter {

    File sumstat_file

    String docker
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String af_col
    String beta_col
    String se_col
    String pval_col

    String outfile = sub(basename(sumstat_file, ".gz"), "\\.bgz$", "") + ".munged.tsv.gz"
    String dollar = "$"

    command <<<

        set -eux

        echo "GWAS meta-analysis - clean and filter sumstats"
        echo "${sumstat_file}"
        echo ""

        catcmd="cat"
        if [[ ${sumstat_file} == *.gz ]] || [[ ${sumstat_file} == *.bgz ]]
        then
            catcmd="zcat"
        fi

        echo "`date` original number of variants"
        $catcmd ${sumstat_file} | tail -n+2 | wc -l

        chr_col=$($catcmd ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "${chr_col}" | head -1 | cut -d ':' -f1)
        pos_col=$($catcmd ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "${pos_col}" | head -1 | cut -d ':' -f1)
        printf "`date` col CHR "${dollar}{chr_col}" col POS "${dollar}{pos_col}"\n"

        $catcmd ${sumstat_file} | awk ' \
            BEGIN{FS="\t| "; OFS="\t"}
            NR==1 {
                for (i=1;i<=NF;i++) {
                    sub("^${chr_col}$", "#CHR", $i);
                    sub("^${pos_col}$", "POS", $i);
                    sub("^${ref_col}$", "REF", $i);
                    sub("^${alt_col}$", "ALT", $i);
                    sub("^${af_col}$", "af_alt", $i);
                    sub("^${beta_col}$", "beta", $i);
                    sub("^${se_col}$", "sebeta", $i);
                    sub("^${pval_col}$", "pval", $i);
                    a[$i]=i;
                    if ($i=="POS") pos=i
                }
                print $0
            } NR>1 {
                sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]); sub("^Y", "24", $a["#CHR"]);
                if ($a["#CHR"] ~ /^[0-9]+$/ && $a["pval"] != 0 && $a["beta"] < 1e6 && $a["beta"] > -1e6 && $a["af_alt"]>0 && (1-$a["af_alt"])>0) {
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
        sort -k$chr_col,${dollar}{chr_col}g -k$pos_col,${dollar}{pos_col}g | \
        bgzip > ${outfile}
        tabix -s $chr_col -b $pos_col -e $pos_col ${outfile}

        echo "`date` new number of variants"
        gunzip -c ${outfile} | tail -n+2 | wc -l
        echo "`date` headers"
        gunzip -c ${outfile} | head -1 | tr '\t' '\n'

        tabix -l ${outfile} > chr.tmp
        echo "`date` $(wc -l chr.tmp | cut -d' ' -f1) chromosomes"
        cat chr.tmp

        echo "`date` unique number of fields"
        gunzip -c ${outfile} | awk 'BEGIN{FS="\t"} {print NF}' | sort -u > n.tmp
        cat n.tmp
        if [ $(wc -l n.tmp | cut -d' ' -f1) != 1 ]; then echo "file not square"; exit 1; fi
        if [ $(wc -l chr.tmp | cut -d' ' -f1) -lt 22 ]; then echo "less than 22 chromosomes"; exit 1; fi

        echo "`date` done"

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# liftover to 38 if needed
task lift {

    File sumstat_file

    String docker
    File b37_ref
    File b38_ref

    File tbi_file = sumstat_file + ".tbi"
    String base = basename(sumstat_file)

    command <<<

        set -eux

        echo "GWAS meta-analysis - lift over sumstats if needed"
        echo "${sumstat_file}"
        echo "${b37_ref}"
        echo "${b38_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${tbi_file} ${base}.tbi

        tabix -R ${b37_ref} ${base} | wc -l > b37.txt
        tabix -R ${b38_ref} ${base} | wc -l > b38.txt

        echo "`date` `cat b37.txt` chr 21 variants build 37"
        echo "`date` `cat b38.txt` chr 21 variants build 38"

        if ((`cat b37.txt` == 0 && `cat b38.txt` == 0)); then
            echo "`date` no chr 21 variants found in either build, quitting"
            exit 1
        fi

        if ((`cat b37.txt` > `cat b38.txt`)); then
            echo "`date` lifting to build 38"
            time lift.py -chr "#CHR" -pos POS -ref REF -alt ALT \
            -chain_file hg19ToHg38.over.chain.gz -tmp_path /cromwell_root/ \
            ${base} > ${base}.lift.out 2> ${base}.lift.err
            gunzip -c ${base}.lifted.gz | \
            cut -f2- | awk '
            BEGIN { FS=OFS="\t" }
            NR==1 {
                for (i=1;i<=NF;i++) {
                    sub("^anew_", "b37_", $i)
                    a[$i]=i;
                }
                print $0
            } NR>1 {
                temp=$a["#CHR"]; $a["#CHR"]=$a["b37_chr"]; $a["b37_chr"]=temp;
                temp=$a["POS"]; $a["POS"]=$a["b37_pos"]; $a["b37_pos"]=temp;
                sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]); sub("^Y", "24", $a["#CHR"]);
                if ($a["#CHR"] ~ /^[0-9]+$/) {
                    print $0
                }
            }' | bgzip > ${base}
        else
            echo "`date` presumably already in build 38"
        fi

    >>>

    output {
        File out = base
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Harmonize variants with gnomad reference
task harmonize {

    File sumstat_file

    String docker
    File gnomad_ref
    String options

    String base = basename(sumstat_file, ".tsv.gz")
    String gnomad_ref_base = basename(gnomad_ref)

    command <<<

        set -eux

        echo "GWAS meta-analysis - harmonize sumstats to reference"
        echo "${sumstat_file}"
        echo "${gnomad_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${gnomad_ref} ${gnomad_ref_base}

        echo "`date` harmonizing stats with gnomAD"
        python3 /META_ANALYSIS/scripts/harmonize.py ${base} ${gnomad_ref_base} ${options} \
        | bgzip > ${base}.${gnomad_ref_base}
        
        tabix -s 1 -b 2 -e 2 ${base}.${gnomad_ref_base}
        echo "`date` done"

    >>>

    output {
        File out = base + "." + gnomad_ref_base
        File out_tbi = base + "." + gnomad_ref_base + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G") + size(gnomad_ref, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Make qc plots
task plot {

    File sumstat_file

    String docker
    Int loglog_ylim

    String base = basename(sumstat_file)

    command <<<

        set -euxo pipefail

        gunzip -c ${sumstat_file} | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) { a[$i]=i; if ($i=="#CHR" || $i=="POS" || $i=="pval" || $i=="af_gnomad" || $i=="af_alt") b[i]=1}}
        {sep=""; for(i=1;i<=NF;i++) if (b[i]==1) { printf sep""$i; sep="\t"} printf "\n"}
        ' | bgzip > ${base} && \

        Rscript - <<EOF
        require(ggplot2)
        require(data.table)
        options(bitmapType='cairo')
        data <- fread("${base}")
        png("${base}_AF.png", width=1000, height=1000, units="px")
        p <- ggplot(data, aes_string(x="af_alt", y="af_gnomad")) +
          geom_point(alpha=0.1) +
          theme_minimal(base_size=18)
        print(p)
        dev.off()
        EOF

        /META_ANALYSIS/scripts/qqplot.R --file ${base} --bp_col "POS" --chrcol "#CHR" --pval_col "pval" --loglog_ylim ${loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}
