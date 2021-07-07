workflow liftover {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (sumstat_file in sumstat_files) {
        call sumstat_to_vcf {
            input: sumstat_file=sumstat_file
        }
        call lift {
            input: sumstat_vcf=sumstat_to_vcf.vcf
        }
    }

    output {
        Array[File] lifted_variants_vcf = lift.lifted_variants
        Array[File] rejected_variants_vcf = lift.rejected_variants
    }
}

task sumstat_to_vcf {

    File sumstat_file

    String docker
    String delim
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String af_col

    String base = basename(sumstat_file, ".gz")

    command <<<

        set -euxo pipefail

        echo "`date` converting sumstat to vcf"

        python3 <<EOF | bgzip > ${base}.vcf.gz

        from datetime import date
        import gzip
        from collections import defaultdict

        sumstat = "${sumstat_file}"
        delim = "${delim}"
        chr_col = "${chr_col}"
        pos_col = "${pos_col}"
        ref_col = "${ref_col}"
        alt_col = "${alt_col}"
        af_col = "${af_col}"

        print('##fileformat=VCFv4.0')

        today = date.today()
        print('##filedate=' + today.strftime('%Y%m%d'))

        print('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">')

        header_line = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        print('\t'.join(header_line))

        with gzip.open(sumstat, 'rt') as f:
            h_idx = {h:i for i,h in enumerate(f.readline().strip().split(delim))}
            h_idx = defaultdict(lambda: 1e9, h_idx)

            for line in f:
                s = line.strip().split(delim)
                s = {i:v for i,v in enumerate(s)}
                s = defaultdict(lambda: '.', s)

                chr = s[h_idx[chr_col]]
                pos = s[h_idx[pos_col]]
                ref = s[h_idx[ref_col]]
                alt = s[h_idx[alt_col]]
                id = ':'.join([chr, pos, ref, alt])
                qual = '.'
                filter = '.'
                info = 'AF=' + s[h_idx[af_col]]

                print('\t'.join([chr, pos, id, ref, alt, qual, filter, info]))

        EOF

        tabix -s 1 -b 2 -e 2 ${base}.vcf.gz

    >>>

    output {
        File vcf = base + ".vcf.gz"
        File vcf_tbi = base + ".vcf.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 3*ceil(size(sumstat_file, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}


task lift {

    File sumstat_vcf

    String docker
    File chainfile
    File b38_assembly_fasta
    File b38_assembly_dict

    String base = basename(sumstat_vcf, ".vcf.gz")

    command <<<

        set -euxo pipefail

        echo "`date` lifting to build 38"
        java -jar /usr/picard/picard.jar LiftoverVcf \
            -I ${sumstat_vcf} \
            -O ${base}.GRCh38.vcf \
            --CHAIN ${chainfile} \
            --REJECT rejected_variants.vcf \
            -R ${b38_assembly_fasta} \
            --MAX_RECORDS_IN_RAM 500000
        
        #bgzip ${base}.GRCh38.vcf > ${base}.GRCh38.vcf.gz
        #bgzip rejected_variants.vcf > rejected_variants.vcf.gz

        #tabix -s 1 -b 2 -e 2 ${base}.GRCh38.vcf.gz
        #tabix -s 1 -b 2 -e 2 rejected_variants.vcf.gz

    >>>

    output {
        File lifted_variants = base + ".GRCh38.vcf"
        #File lifted_variants_tbi = base + ".GRCh38.vcf.gz.tbi"
        File rejected_variants = "rejected_variants.vcf"
        #File rejected_variants_tbi = "rejected_variants.vcf.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "32 GB"
        disks: "local-disk " + 10*ceil(size(sumstat_vcf, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}