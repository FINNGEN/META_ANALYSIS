version 1.0

workflow liftover {

    input {
        File sumstats_loc
        Array[String] sumstat_files = read_lines(sumstats_loc)

        String delim
        String chr_col
        String pos_col
        String ref_col
        String alt_col
        String af_col
        String beta_col
    }

    scatter (sumstat_file in sumstat_files) {
        call sumstat_to_vcf {
            input:
                sumstat_file = sumstat_file,
                delim = delim,
                chr_col = chr_col,
                pos_col = pos_col,
                ref_col = ref_col,
                alt_col = alt_col
        }
        call lift {
            input:
                sumstat_vcf = sumstat_to_vcf.vcf
        }
        call lift_postprocess {
            input:
                lifted_vcf = lift.lifted_variants_vcf,
                sumstat_file = sumstat_file,
                delim = delim,
                chr_col = chr_col,
                pos_col = pos_col,
                ref_col = ref_col,
                alt_col = alt_col,
                af_col = af_col,
                beta_col = beta_col
        }
    }

    output {
        Array[File] lifted_variants = lift_postprocess.lifted_variants
        Array[File] rejected_variants_vcf = lift.rejected_variants
    }
}

task sumstat_to_vcf {

    input {
        File sumstat_file
        String delim
        String chr_col
        String pos_col
        String ref_col
        String alt_col

        String docker

        String base = basename(sumstat_file, ".gz")
    }

    command <<<

        set -euxo pipefail

        echo "`date` converting sumstat to vcf"

        python3 <<EOF | bgzip > ~{base}.vcf.gz

        from datetime import date
        import gzip
        from collections import defaultdict

        sumstat = "~{sumstat_file}"
        delim = "~{delim}"
        chr_col = "~{chr_col}"
        pos_col = "~{pos_col}"
        ref_col = "~{ref_col}"
        alt_col = "~{alt_col}"

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

        String base = basename(sumstat_vcf, ".vcf.gz")
    }

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
        String delim
        String chr_col
        String pos_col
        String ref_col
        String alt_col
        String af_col
        String beta_col

        String docker

        String base = basename(lifted_vcf, ".vcf")
    }

    command <<<

        set -euxo pipefail

        # Sort by old positions
        grep -v "^#" ~{lifted_vcf} | tr ":" "\t" | awk '
        BEGIN{OFS="\t"}
        { gsub("chr", ""); gsub("X", "23"); gsub("Y", "24"); print $1,$2,$7,$8,$3,$4,$5,$6,$11 }
        ' | sort -k5,5g -k6,6g > ~{base}.tsv

        python3 <<EOF | sort -k1,1g -k2,2g | bgzip > ~{base}.tsv.gz

        import gzip
        from collections import defaultdict

        valid_chrs = [str(i) for i in range(1,25)]

        sumstat = "~{sumstat_file}"
        delim = "~{delim}"
        chr_col = "~{chr_col}"
        pos_col = "~{pos_col}"
        ref_col = "~{ref_col}"
        alt_col = "~{alt_col}"
        af_col = "~{af_col}"
        beta_col = "~{beta_col}"

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
                    except ValueError:
                        break
                    sumstat_ref = sumstat_line[sumstat_h_idx[ref_col]]
                    sumstat_alt = sumstat_line[sumstat_h_idx[alt_col]]

        EOF

        tabix -S 1 -s 1 -b 2 -e 2 ~{base}.tsv.gz

    >>>

    output {
        File lifted_variants = base + ".tsv.gz"
        File lifted_variants_tbi = base + ".tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 4*ceil(size(lifted_vcf, "G") + size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}