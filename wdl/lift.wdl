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
        String? af_col
        String? beta_col
    }

    scatter (sumstat_file in sumstat_files) {
        call sumstat_to_vcf {
            input:
                sumstat_file = sumstat_file,
                delim = delim,
                chr_col = chr_col,
                pos_col = pos_col,
                ref_col = ref_col,
                alt_col = alt_col,
                af_col = af_col,
                beta_col = beta_col
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
        String? af_col
        String? beta_col

        String docker

        String base = basename(sumstat_file, ".gz")
    }

    command <<<

        set -euxo pipefail

        echo "`date` converting sumstat to vcf"

        python3 <<EOF | bgzip > ~{base}.vcf.gz

        from datetime import date
        import gzip
        import sys
        from collections import defaultdict
        from contextlib import contextmanager

        @contextmanager
        def uopen(fname,oper_type,encoding="utf-8"):
            """
            Universal opener
            Open both gzipped and plaintext files with no fuzz
            """
            gz_magicnumber=b"\x1f\x8b"
            type="normal"
            with open(fname,"rb") as f:
                if f.read(2) == gz_magicnumber:
                    type="gz"
            if type=="normal":
                with open(fname,oper_type,encoding=encoding) as f:
                    yield f
            elif type == "gz":
                with gzip.open(fname,oper_type,encoding=encoding) as f:
                    yield f
            else:
                raise Exception("invalid file format")

        sumstat = "~{sumstat_file}"
        delim = "~{delim}"
        chr_col = "~{chr_col}".lstrip('#')
        pos_col = "~{pos_col}".lstrip('#')
        ref_col = "~{ref_col}"
        alt_col = "~{alt_col}"
        af_col_str = ~{if defined(af_col) then "'~{af_col}'" else "None"}
        beta_col_str = ~{if defined(beta_col) then "'~{beta_col}'" else "None"}

        # Parse comma-separated lists of columns for allele frequency and beta values, if provided
        af_cols = [col.strip() for col in af_col_str.split(',')] if af_col_str is not None else []
        beta_cols = [col.strip() for col in beta_col_str.split(',')] if beta_col_str is not None else []

        with uopen(sumstat, 'rt') as f:
            sumstat_header = f.readline().strip().split(delim)
            # Remove leading '#' from first column if present, for consistent lookups
            if sumstat_header and sumstat_header[0].startswith('#'):
                sumstat_header[0] = sumstat_header[0].lstrip('#')
            header_set = set(sumstat_header)

            # Validate that every required input column exists in the sumstat header
            required_cols = [chr_col, pos_col, ref_col, alt_col] + af_cols + beta_cols
            missing_cols = [c for c in required_cols if c not in header_set]
            if missing_cols:
                sys.stderr.write(
                    'ERROR: the following input columns were not found in the sumstat header of ' + sumstat + ':\n'
                    + '  missing: ' + ', '.join(missing_cols) + '\n'
                    + '  available: ' + ', '.join(sumstat_header) + '\n'
                )
                sys.exit(1)

            h_idx = {h:i for i,h in enumerate(sumstat_header)}
            h_idx = defaultdict(lambda: 1e9, h_idx)

        print('##fileformat=VCFv4.0')

        today = date.today()
        print('##filedate=' + today.strftime('%Y%m%d'))

        header_line = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        print('\t'.join(header_line))

        with uopen(sumstat, 'rt') as f:
            f.readline()  # skip header (already validated above)

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
            --RECOVER_SWAPPED_REF_ALT true \
            --DISABLE_SORT true

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
        String? af_col
        String? beta_col

        String docker

        String base = basename(lifted_vcf, ".vcf")
    }

    command <<<

        set -eux

        # Sort by old positions
        grep -v "^#" ~{lifted_vcf} | tr ":" "\t" | awk '
        BEGIN{OFS="\t"}
        { gsub("chr", ""); gsub("X", "23"); gsub("Y", "24"); print $1,$2,$7,$8,$3,$4,$5,$6,$11 }
        ' | sort -t'~{delim}' -k5,5g -k6,6g > ~{base}.tsv

        python3 <<EOF | bgzip > ~{base}.unsorted.tsv.gz

        import gzip
        from collections import defaultdict
        from contextlib import contextmanager

        @contextmanager
        def uopen(fname,oper_type,encoding="utf-8"):
            """
            Universal opener
            Open both gzipped and plaintext files with no fuzz
            """
            gz_magicnumber=b"\x1f\x8b"
            type="normal"
            with open(fname,"rb") as f:
                if f.read(2) == gz_magicnumber:
                    type="gz"
            if type=="normal":
                with open(fname,oper_type,encoding=encoding) as f:
                    yield f
            elif type == "gz":
                with gzip.open(fname,oper_type,encoding=encoding) as f:
                    yield f
            else:
                raise Exception("invalid file format")

        valid_chrs = [str(i) for i in range(1,25)]

        sumstat = "~{sumstat_file}"
        delim = "~{delim}"
        chr_col = "~{chr_col}".lstrip('#')
        pos_col = "~{pos_col}".lstrip('#')
        ref_col = "~{ref_col}"
        alt_col = "~{alt_col}"
        af_col_str = ~{if defined(af_col) then "'~{af_col}'" else "None"}
        beta_col_str = ~{if defined(beta_col) then "'~{beta_col}'" else "None"}

        # Parse comma-separated lists of columns for allele frequency and beta values, if provided
        af_cols = [col.strip() for col in af_col_str.split(',')] if af_col_str is not None else []
        beta_cols = [col.strip() for col in beta_col_str.split(',')] if beta_col_str is not None else []

        with uopen(sumstat, 'rt') as s_f:
            sumstat_header = s_f.readline().strip().split(delim)
            # Remove leading '#' from first column if present, for consistent indexing
            if sumstat_header and sumstat_header[0].startswith('#'):
                sumstat_header[0] = sumstat_header[0].lstrip('#')
            sumstat_h_idx = {h:i for i,h in enumerate(sumstat_header)}
            sumstat_h_idx = defaultdict(lambda: None, sumstat_h_idx)
            sumstat_header.extend(['b37_chr', 'b37_pos', 'b37_ref', 'b37_alt', 'liftover_info'])
            # Ensure header starts with '#' for tabix compatibility
            header_line = delim.join(sumstat_header)
            if not header_line.startswith('#'):
                header_line = '#' + header_line
            print(header_line)
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
                            for af_col in af_cols:
                                sumstat_af = sumstat_line[sumstat_h_idx[af_col]] if sumstat_line[sumstat_h_idx[af_col]] != "NA" else None
                                sumstat_line[sumstat_h_idx[af_col]] = str(1 - float(sumstat_af) if sumstat_af is not None else "NA")
                            for beta_col in beta_cols:
                                sumstat_beta = sumstat_line[sumstat_h_idx[beta_col]] if sumstat_line[sumstat_h_idx[beta_col]] != "NA" else None
                                sumstat_line[sumstat_h_idx[beta_col]] = str(-1 * float(sumstat_beta) if sumstat_beta is not None else "NA")
                        sumstat_line.extend([str(old_chr), str(old_pos), old_ref, old_alt, info])
                        print(delim.join(sumstat_line))
                        sumstat_line = s_f.readline().strip().split(delim)
                        if len(sumstat_line) == 1: # final newline
                            break
                        try:
                            sumstat_chr = int(sumstat_line[sumstat_h_idx[chr_col]].replace('chr', '').replace('X', '23').replace('Y', '24'))
                            sumstat_pos = int(sumstat_line[sumstat_h_idx[pos_col]])
                        except ValueError:
                            break
                        sumstat_ref = sumstat_line[sumstat_h_idx[ref_col]]
                        sumstat_alt = sumstat_line[sumstat_h_idx[alt_col]]

        EOF

        # Strip leading '#' from column names for searching
        chr_col_clean=`echo "~{chr_col}" | sed 's/^#//'`
        pos_col_clean=`echo "~{pos_col}" | sed 's/^#//'`
        
        chr_col_n=`zcat ~{base}.unsorted.tsv.gz | head -n 1 | sed 's/^#//' | awk -v RS='~{delim}' -v col="$chr_col_clean" '$0 == col {print NR; exit}'`
        pos_col_n=`zcat ~{base}.unsorted.tsv.gz | head -n 1 | sed 's/^#//' | awk -v RS='~{delim}' -v col="$pos_col_clean" '$0 == col {print NR; exit}'`

        cat \
        <(zcat ~{base}.unsorted.tsv.gz | sed '1!d') \
        <(zcat ~{base}.unsorted.tsv.gz | tail -n+2 | sort -t'~{delim}' -k${chr_col_n},${chr_col_n}V -k${pos_col_n},${pos_col_n}g) | \
        bgzip > ~{base}.tsv.gz

        if [[ "~{delim}" == $'\t' ]]; then
          tabix -s ${chr_col_n} -b ${pos_col_n} -e ${pos_col_n} ~{base}.tsv.gz
        else
          touch ~{base}.tsv.gz.tbi
        fi

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
