task IgCaller {

    String id
    File tumor_bam
    File tumor_bai

    #Need normal or reference
    File? normal_bam
    File? normal_bai
    File? reference_genome

    #Other IgCaller args
    String? igc_ref_file_path="/src/IgCaller/IgCaller_reference_files/"
    String? assembly="hg19"
    String? chrAnnotation="ensembl"
    String? keep_minibam="no"
    String? seqtype="wgs"
    String? other_args=""

    #Runtime args
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/igcaller:1.1"
    Int? disk_size = ceil(size(tumor_bam, "GB") + size(normal_bam, "GB"))
    Int? disk_size_add ="50" #ceil(size(tumor_bam, "GB") / 2)
    Int? disk_size_total = disk_size + disk_size_add
    Int? threads="1"
    Float? memory="3.75"

    command {
        set -euo pipefail

        #Change file names to get output files in wanted format
        mkdir input_bams
        mv ${tumor_bam} input_bams/${id}.bam
        mv ${tumor_bai} input_bams/${id}.bam.bai

        ### Prep args (need to do here, because STAR output is piped into Arriba)
        all_other_args="${other_args}"

        if [ -n "${reference_genome}" ]; then
            all_other_args="$all_other_args -R ${reference_genome}"
        fi

        if [ -n "${normal_bam}" ]; then
            mv ${normal_bam} input_bams/${id}_normal.bam
            mv ${normal_bai} input_bams/${id}_normal.bam.bai
            all_other_args="$all_other_args -N input_bams/${id}_normal.bam"
        fi

        mkdir igcaller_outdir

        IgCaller -I ${igc_ref_file_path} \
            -V ${assembly} \
            -C ${chrAnnotation} \
            -T input_bams/${id}.bam \
            -o igcaller_outdir \
            -kmb ${keep_minibam} \
            -seq ${seqtype} \
            $all_other_args

            # Move output for convenience
            # gzip files that are less likely to be used
            mv igcaller_outdir/*/* .
            gzip ${id}_output_IGH.tsv
            gzip ${id}_output_IGK.tsv
            gzip ${id}_output_IGL.tsv
            gzip ${id}_output_CSR.tsv
    }

    output {
        File igcaller_filtered_tsv="${id}_output_filtered.tsv"
        File igcaller_IGH="${id}_output_IGH.tsv.gz"
        File igcaller_IGK="${id}_output_IGK.tsv.gz"
        File igcaller_IGL="${id}_output_IGL.tsv.gz"
        File igcaller_CSR="${id}_output_CSR.tsv.gz"
        File igcaller_IG_rearrangements="${id}_output_oncogenic_IG_rearrangements.tsv"

    }

    runtime {
        docker: "${docker}"
        disks: "local-disk ${disk_size_total} HDD"
        memory: "${memory}GB"
        cpu: "${threads}"
        preemptible: 4
    }

    meta {
        author : "bknisbac"
        email : "bknisbac@broadinstitute.org"
    }
}

workflow IgCaller_workflow {
    call IgCaller
}
