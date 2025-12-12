version 1.0

workflow munge {

    input {
        File sumstats_loc
        Boolean liftover
        String common_docker
        String meta_docker
    }

    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        
        call clean_filter {
            input:
            sumstat_file = sumstat_file,
            docker = meta_docker
        }

        if (liftover) {
            call sumstat_to_vcf {
                input:
                    sumstat_file = clean_filter.out,
                    docker = common_docker
            }
            
            call lift {
                input:
                    sumstat_vcf = sumstat_to_vcf.vcf
            }
            
            call lift_postprocess {
                input:
                    lifted_vcf = lift.lifted_variants_vcf,
                    sumstat_file = clean_filter.out,
                    docker = common_docker
            }
        }
        
        call harmonize {
            input:
                sumstat_file = select_first([lift_postprocess.lifted_variants, clean_filter.out]),
                docker = meta_docker
        }
        
        call filt {
            input:
                sumstat_file = harmonize.out,
                docker = meta_docker
        }
    }

    output {
        Array[File] cleaned_sumstats = clean_filter.out
        Array[File] harmonized_sumstats = harmonize.out
        Array[File] filtered_sumstats = filt.out
        Array[Array[File]] plots = filt.pngs
    }
}

# Filter bad quality variants
task clean_filter {

    input {
        File sumstat_file

        String docker
        String chr_col
        String pos_col
        String ref_col
        String alt_col
        String af_col
        String effect_col
        String se_col
        String pval_col
        String options
    }

    String outfile = sub(basename(sumstat_file, ".gz"), "\\.bgz$", "") + ".munged.tsv.gz"

    command <<<

        set -eux

        echo "GWAS meta-analysis - clean and filter sumstats"
        echo "~{sumstat_file}"
        echo ""

        echo "`date` original number of variants"
        zcat -f ~{sumstat_file} | tail -n+2 | wc -l

        chr_col=$(zcat -f ~{sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "~{chr_col}" | head -1 | cut -d ':' -f1)
        pos_col=$(zcat -f ~{sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "~{pos_col}" | head -1 | cut -d ':' -f1)
        printf "`date` col CHR "${chr_col}" col POS "${pos_col}"\n"

        munge.py ~{sumstat_file} \
            --chr-col "~{chr_col}" \
            --pos-col "~{pos_col}" \
            --ref-col "~{ref_col}" \
            --alt-col "~{alt_col}" \
            --af-col "~{af_col}" \
            --effect-col "~{effect_col}" \
            --se-col "~{se_col}" \
            --pval-col "~{pval_col}" \
            ~{options} | \
            sort -k$chr_col,${chr_col}g -k$pos_col,${pos_col}g | \
            bgzip > ~{outfile}
        
        tabix -S 1 -s $chr_col -b $pos_col -e $pos_col ~{outfile}

        echo "`date` headers"
        gunzip -c ~{outfile} | head -1 | tr '\t' '\n' | nl

        tabix -l ~{outfile} > chr.tmp
        echo "`date` $(wc -l chr.tmp | cut -d' ' -f1) chromosomes"
        cat chr.tmp

        if [ $(wc -l chr.tmp | cut -d' ' -f1) -lt 22 ]; then echo "less than 22 chromosomes"; exit 1; fi

        echo "`date` done"

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Convert sumstat to vcf for liftover
task sumstat_to_vcf {

    input {
        File sumstat_file

        String docker
    }

    String base = basename(sumstat_file, ".gz")

    command <<<

        set -euxo pipefail

        echo "`date` converting sumstat to vcf"

        python3 <<EOF | bgzip > ~{base}.vcf.gz

        from datetime import date
        import gzip
        from collections import defaultdict

        sumstat = '~{sumstat_file}'
        delim = '\t'
        chr_col = '#CHR'
        pos_col = 'POS'
        ref_col = 'REF'
        alt_col = 'ALT'

        print('##fileformat=VCFv4.0')

        today = date.today()
        print('##filedate=' + today.strftime('%Y%m%d'))

        header_line = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        print('\t'.join(header_line))

        with gzip.open(sumstat, 'rt') as f:
            h_idx = {h:i for i,h in enumerate(f.readline().strip().split(delim))}
            h_idx = defaultdict(lambda: 1e9, h_idx)

            for line in f:
                s = line.strip().split(delim)
                s = {i:v for i,v in enumerate(s)}
                s = defaultdict(lambda: '.', s)

                chr = str(s[h_idx[chr_col]])
                if chr == '23':
                    chr = 'X'
                if chr == '24':
                    chr = 'Y'
                if chr == '25':
                    chr = 'M'
                if chr[:3] != 'chr':
                    chr = 'chr' + chr

                pos = s[h_idx[pos_col]]
                ref = s[h_idx[ref_col]]
                alt = s[h_idx[alt_col]]
                id = ':'.join([chr, pos, ref, alt])
                qual = '.'
                filter = '.'
                info = '.'

                print('\t'.join([chr, pos, id, ref, alt, qual, filter, info]))

        EOF

        tabix -s 1 -b 2 -e 2 ~{base}.vcf.gz

    >>>

    output {
        File vcf = base + ".vcf.gz"
        File vcf_tbi = base + ".vcf.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 3*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}


task lift {

    input {
        File sumstat_vcf

        String docker
        File chainfile
        File b38_assembly_fasta
        File b38_assembly_dict
    }

    String base = basename(sumstat_vcf, ".vcf.gz")

    command <<<

        set -euxo pipefail

        echo "`date` lifting to build 38"
        java -jar /usr/picard/picard.jar LiftoverVcf \
            -I ~{sumstat_vcf} \
            -O ~{base}.GRCh38.vcf \
            --CHAIN ~{chainfile} \
            --REJECT rejected_variants.vcf \
            -R ~{b38_assembly_fasta} \
            --MAX_RECORDS_IN_RAM 500000 \
            --RECOVER_SWAPPED_REF_ALT true

    >>>

    output {
        File lifted_variants_vcf = base + ".GRCh38.vcf"
        File rejected_variants = "rejected_variants.vcf"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "32 GB"
        disks: "local-disk " + 10*ceil(size(sumstat_vcf, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}


task lift_postprocess {

    input {
        File lifted_vcf
        File sumstat_file

        String docker
    }

    String base = basename(lifted_vcf, ".vcf")

    command <<<

        set -eux

        # Sort by old positions
        grep -v "^#" ~{lifted_vcf} | tr ":" "\t" | awk '
        BEGIN{OFS="\t"}
        { gsub("chr", ""); gsub("X", "23"); gsub("Y", "24"); print $1,$2,$7,$8,$3,$4,$5,$6,$11 }
        ' | sort -k5,5g -k6,6g > ~{base}.tsv

        chr_col=$(zcat ~{sumstat_file} | head -1 | tr '\t ' '\n' | grep -nwF "#CHR" | head -1 | cut -d ':' -f1)
        pos_col=$(zcat ~{sumstat_file} | head -1 | tr '\t ' '\n' | grep -nwF "POS" | head -1 | cut -d ':' -f1)

        python3 <<EOF | sort -k$chr_col,${chr_col}g -k$pos_col,${pos_col}g | bgzip > ~{base}.tsv.gz

        import gzip
        from collections import defaultdict

        valid_chrs = set([str(i) for i in range(1,25)])

        sumstat = '~{sumstat_file}'
        delim = '\t'
        chr_col = '#CHR'
        pos_col = 'POS'
        ref_col = 'REF'
        alt_col = 'ALT'
        af_col = 'af_alt'
        beta_col = 'beta'

        s_f = gzip.open(sumstat, 'rt')
        sumstat_header = s_f.readline().strip().split(delim)
        sumstat_h_idx = {h:i for i,h in enumerate(sumstat_header)}
        sumstat_h_idx = defaultdict(lambda: None, sumstat_h_idx)
        sumstat_header.extend(['b37_chr', 'b37_pos', 'b37_ref', 'b37_alt', 'liftover_info'])
        print(delim.join(sumstat_header))
        sumstat_chr = 0
        sumstat_pos = 0
        sumstat_ref = 'N'
        sumstat_alt = 'N'

        with open('~{base}.tsv', 'rt') as f:
            for line in f:
                s = line.strip().split('\t')
                new_chr = s[0]
                if new_chr not in valid_chrs:
                    continue
                new_pos = s[1]
                new_ref = s[2]
                new_alt = s[3]
                old_chr = int(s[4])
                old_pos = int(s[5])
                old_ref = s[6]
                old_alt = s[7]
                info = s[8]
                info_list = info.strip().split(';')

                while sumstat_chr < old_chr or (sumstat_chr == old_chr and sumstat_pos < old_pos):
                    sumstat_line = s_f.readline().strip().split(delim)
                    sumstat_chr = int(sumstat_line[sumstat_h_idx[chr_col]].replace('chr', '').replace('X', '23').replace('Y', '24'))
                    sumstat_pos = int(sumstat_line[sumstat_h_idx[pos_col]])
                    sumstat_ref = sumstat_line[sumstat_h_idx[ref_col]]
                    sumstat_alt = sumstat_line[sumstat_h_idx[alt_col]]
                if sumstat_chr == old_chr and sumstat_pos == old_pos and sumstat_ref == old_ref and sumstat_alt == old_alt:
                    sumstat_line[sumstat_h_idx[chr_col]] = new_chr
                    sumstat_line[sumstat_h_idx[pos_col]] = new_pos
                    sumstat_line[sumstat_h_idx[ref_col]] = new_ref
                    sumstat_line[sumstat_h_idx[alt_col]] = new_alt
                    if 'SwappedAlleles' in info_list:
                        sumstat_af = sumstat_line[sumstat_h_idx[af_col]] if sumstat_line[sumstat_h_idx[af_col]] != "NA" else None
                        sumstat_beta = sumstat_line[sumstat_h_idx[beta_col]] if sumstat_line[sumstat_h_idx[beta_col]] != "NA" else None
                        sumstat_line[sumstat_h_idx[af_col]] = str(1 - float(sumstat_af) if sumstat_af is not None else "NA")
                        sumstat_line[sumstat_h_idx[beta_col]] = str(-1 * float(sumstat_beta) if sumstat_beta is not None else "NA")
                    sumstat_line.extend([str(old_chr), str(old_pos), old_ref, old_alt, info])
                    print(delim.join(sumstat_line))
                    sumstat_line = s_f.readline().strip().split(delim)
                    try:
                        sumstat_chr = int(sumstat_line[sumstat_h_idx[chr_col]].replace('chr', '').replace('X', '23').replace('Y', '24'))
                        sumstat_pos = int(sumstat_line[sumstat_h_idx[pos_col]])
                    except (IndexError, ValueError):
                        break
                    sumstat_ref = sumstat_line[sumstat_h_idx[ref_col]]
                    sumstat_alt = sumstat_line[sumstat_h_idx[alt_col]]

        EOF

        tabix -s $chr_col -b $pos_col -e $pos_col ~{base}.tsv.gz

    >>>

    output {
        File lifted_variants = base + ".tsv.gz"
        File lifted_variants_tbi = base + ".tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + (5 * (ceil(size(lifted_vcf, "G") + size(sumstat_file, "G")) + 1)) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Harmonize variants with gnomad reference
task harmonize {

    input {
        File sumstat_file

        String docker
        File gnomad_ref
        String options
    }

    String base = basename(sumstat_file, ".tsv.gz")
    String gnomad_ref_base = basename(gnomad_ref)

    command <<<

        set -eux

        echo "GWAS meta-analysis - harmonize sumstats to reference"
        echo "~{sumstat_file}"
        echo "~{gnomad_ref}"
        echo ""

        mv ~{sumstat_file} ~{base}
        mv ~{gnomad_ref} ~{gnomad_ref_base}

        echo "`date` harmonizing stats with gnomAD"
        harmonize.py ~{base} ~{gnomad_ref_base} ~{options} | bgzip > ~{base}.~{gnomad_ref_base}
        
        tabix -s 1 -b 2 -e 2 ~{base}.~{gnomad_ref_base}
        echo "`date` done"

    >>>

    output {
        File out = base + "." + gnomad_ref_base
        File out_tbi = base + "." + gnomad_ref_base + ".tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G") + size(gnomad_ref, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

# Make qc plots
task filt {

    input {
        File sumstat_file

        String docker
        Int loglog_ylim
        Float max_af_diff
    }

    String base = basename(sumstat_file)

    command <<<

        set -euxo pipefail

        mv ~{sumstat_file} ~{base}

        Rscript - <<EOF
        require(ggplot2)
        require(data.table)
        options(bitmapType='cairo')
        data <- fread("~{base}", select=c("af_gnomad", "af_alt"), header=T)
        message(nrow(data[abs(af_alt - af_gnomad) > ~{max_af_diff}]), " variants will be filtered out with AF difference to GnomAD > ~{max_af_diff}")
        png("~{base}_AF.png", width=1000, height=1000, units="px")
        p <- ggplot(data, aes(x = af_alt, y = af_gnomad, color = abs(af_alt - af_gnomad) > .2)) +
          geom_point(alpha=0.1) +
          theme_minimal(base_size=18) +
          scale_color_manual(values = c("black", "red")) +
          theme(legend.position = "none")
        print(p)
        dev.off()
        EOF

        zcat -f ~{base} | awk '
        function abs(x) { return ((x < 0.0) ? -x : x) }
        BEGIN { FS=OFS="\t" }
        NR==1 { for(i=1;i<=NF;i++) h[$i]=i; print }
        NR>1 && (abs($h["af_alt"] - $h["af_gnomad"]) <= ~{max_af_diff} || $h["af_gnomad"] == "NA") { print }
        ' | bgzip > tmp

        mv tmp ~{base}
        tabix -s 1 -b 2 -e 2 ~{base}

        qqplot.R --file ~{base} --bp_col "POS" --chrcol "#CHR" --pval_col "pval" --loglog_ylim ~{loglog_ylim}

    >>>

    output {
        Array[File] pngs = glob("*.png")
        File out = base
        File out_tbi = base + ".tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "20 GB"
        disks: "local-disk " + 5*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}
