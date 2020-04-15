
task liftfile {
    File f
    File chainfile
    String chr_col
    String pos_col
    String ref_col
    String alt_col

    String docker

    String base=basename(f)

    command <<<
        lift.py -chr ${chr_col} -pos ${pos_col} -ref ${ref_col} -alt ${alt_col} \
        -chain_file ${chainfile} ${f} -tmp_path /cromwell_root/
    >>>
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }

    output {
        File lifted="${base}.lifted.gz"
        File lifted_idx="${base}.lifted.gz.tbi"
        Array[File] lifterrors=glob("*_errors")
    }
}

workflow lift {
    File files
    Array[String] fs= read_lines(files)
    scatter (f in fs) {
        call liftfile {
            input:  f=f
        }
    }

}
