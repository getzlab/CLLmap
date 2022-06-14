
task minibam {

    String sample_id
    File bam
    File? rule_file
    String? regions_opt = "-l" #-l gets mate too (see specs for other options: -g, -G, -L, -k)
    File? regions = "gs://bknisbac1/resources/ig_regions.hg19.bed"
	String? regions_str_opt = "-l"
    String? regions_str = ""
    String? command_line = ""
    String? minibam_descriptor = "mini"
    Int? region_pad = "0"

    Int? memory = "15"
    Int? disk_space_add = "50"
    Int? disk_space = ceil(size(bam ,"GB")) + disk_space_add
    Int? num_threads = "4"
    Int? num_preempt = "4"

    Int? num_threads_sort = num_threads - 1
    Int? memory_per_thread_sort = memory / num_threads

    command {
        set -euo pipefail

        #Define output prefix
        minibam_prefix=${sample_id}.${minibam_descriptor}

        # placeholders for optional outputs
        touch $minibam_prefix + ".bam.bai"
        #touch $minibam_prefix + ".bam.idxstats"

        ### VariantBam command and parametres ###
        variant_args=" -v -b -t ${num_threads}"
        #Rule file for selection
        if [ -n "${rule_file}" ]; then
            variant_args = "$variant_args --rules ${rule_file}"
        fi
        #Regions to analyze (BED format)
        if [ -n "${regions}" ]; then
            variant_args="$variant_args ${regions_opt} ${regions}"
        fi
        #Regions to analyze (samtools format string, e.g. 14:1323-1329)
        if [ -n "${regions_str}" ]; then
            variant_args="$variant_args ${regions_str_opt} ${regions_str}"
        fi
        #Pad regions
        if [ -n "${region_pad}" ]; then
            variant_args="$variant_args -P ${region_pad}"
        fi
        #Free-style command line addition
        if [ -n "${command_line}" ]; then
            variant_args="$variant_args ${command_line}"
        fi
        #Run command
        variant ${bam} $variant_args > $minibam_prefix.bam

        ### Not needed in this specific workflow
        ### sort + index (optional) ###
        #samtools sort -O BAM -@ "${num_threads_sort}" -m ${memory_per_thread_sort}"G" -o $minibam_prefix.bam $minibam_prefix.bam
        #samtools index $minibam_prefix.bam

        ### idxstats ###
        #samtools idxstats $minibam_prefix.bam > $minibam_prefix.bam.idxstats
    }

    output {
        File minibam = "${sample_id}.${minibam_descriptor}.bam"
    }

    runtime {
        docker: "gcr.io/broad-cga-bknisbac-wupo1/variantbam:1.0"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Binyamin A. Knisbacher"
    }
}


task samtofastq {

    File input_bam_cram
    String prefix
    File? reference_fa

    Float? memory = "3.75"
    Int java_memory = floor(memory - 0.5)
    Int? disk_space = ceil(size(input_bam_cram, "G")*2 + size(reference_fa, "G"))
    Int? num_threads = "1"
    Int? num_preempt = "4"

    command {
        set -euo pipefail

        # make sure path is absolute
        input_bam_abs=${input_bam_cram}
        if [[ $input_bam_abs != /* ]]; then
            input_bam_abs=$PWD/$input_bam_abs
        fi

        mkdir samtofastq  # workaround for named pipes
        python3 -u /src/run_SamToFastq.py $input_bam_abs -p ${prefix} ${"--reference_fa " + reference_fa} --output_dir samtofastq --memory ${java_memory}
        mv samtofastq/${prefix}_*.fastq.gz .
    }

    output {
        File fastq1="${prefix}_1.fastq.gz"
        File fastq2="${prefix}_2.fastq.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}

task bwa_mem {
	File fastq1
    File fastq2
    String sample_id
    File reference_fa
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_dict
    File reference_idx
    File reference_pac
    File reference_sa

    Int threads_bwa="8"
    Int threads_sort="4"
    Int num_preempt="4"
    Int mem="16"
    Int disk_size = ceil(size(fastq1, "G") + size(fastq2, "G") + size(reference_fa, "G")*2) * 2
    Int threads = threads_bwa + threads_sort
    command <<<
      set -euo pipefail
      /bwa/bwa mem -t ${threads_bwa} ${reference_fa} ${fastq1} ${fastq2} -R "@RG\tID:${sample_id}\tSM:${sample_id}_${sample_id}\tLB:${sample_id}_${sample_id}\tPL:ILLUMINA" | samtools sort -@${threads_sort} -o ${sample_id}.bam -
      samtools index -@ ${threads} ${sample_id}.bam
      samtools idxstats ${sample_id}.bam > ${sample_id}.bam.idxstats.txt
    >>>

    output{
    	File bam_out = "${sample_id}.bam"
        File bam_index = "${sample_id}.bam.bai"
        File bam_idxstats = "${sample_id}.bam.idxstats.txt"
    }

    runtime {
        docker: "alexbarbera/samtools_bwa:v1.0"
        memory: "${mem}GB"
        disks: "local-disk ${disk_size} HDD"
        cpu: "${threads}"
        preemptible: "${num_preempt}"
    }

    meta{
        #Adapted partially from Alex Barbera
        author: "Binyamin Knisbacher"
    }
}

task IgCaller {

    String id
    File tumor_bam
    File tumor_bai

    #Need normal or reference
    File? normal_bam
    File? normal_bai
    File? reference_fa

    #Other IgCaller args
    String? igc_ref_file_path="/src/IgCaller/IgCaller_reference_files/"
    String? assembly="hg19"
    String? chrAnnotation="ensembl"
    String? keep_minibam="no"
    String? seqtype="wgs"
    String? other_args=""

    #Runtime args
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/igcaller:1.1c"
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

        if [ -n "${reference_fa}" ]; then
            all_other_args="$all_other_args -R ${reference_fa}"
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
        author : "Binyamin Knisbacher"
        email : "bknisbac@broadinstitute.org"
    }
}


workflow igcaller_realign_workflow {

    String pair_id
    String tumor_sample_id
    String normal_sample_id
    File tumor_bam
    File normal_bam
    File reference_fa

    #maybe a few unneeded indexes here (only needed for bwa-mem)
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_dict
    File reference_idx
    File reference_pac
    File reference_sa

    String? minibam_descriptor = 'minibam_IG'

    call minibam as minibam_tumor {
        input:
            sample_id=tumor_sample_id,
            bam=tumor_bam,
            minibam_descriptor=minibam_descriptor
    }

    call minibam as minibam_normal {
        input:
            sample_id=normal_sample_id,
            bam=normal_bam,
            minibam_descriptor=minibam_descriptor
    }

    call samtofastq as samtofastq_tumor {
        input:
            prefix=tumor_sample_id,
            input_bam_cram=minibam_tumor.minibam,
            reference_fa=reference_fa
    }

    call samtofastq as samtofastq_normal {
        input:
            prefix=normal_sample_id,
            input_bam_cram=minibam_normal.minibam,
            reference_fa=reference_fa
    }

    call bwa_mem as realign_tumor {
        input:
            fastq1=samtofastq_tumor.fastq1,
            fastq2=samtofastq_tumor.fastq2,
            sample_id=tumor_sample_id,
            reference_fa=reference_fa,
            reference_amb=reference_amb,
            reference_ann=reference_ann,
            reference_bwt=reference_bwt,
            reference_dict=reference_dict,
            reference_idx=reference_idx,
            reference_pac=reference_pac,
            reference_sa=reference_sa
    }

    call bwa_mem as realign_normal {
        input:
            fastq1=samtofastq_normal.fastq1,
            fastq2=samtofastq_normal.fastq2,
            sample_id=normal_sample_id,
            reference_fa=reference_fa,
            reference_amb=reference_amb,
            reference_ann=reference_ann,
            reference_bwt=reference_bwt,
            reference_dict=reference_dict,
            reference_idx=reference_idx,
            reference_pac=reference_pac,
            reference_sa=reference_sa
    }

    call IgCaller {
        input:
            id=pair_id,
            tumor_bam=realign_tumor.bam_out,
            tumor_bai=realign_tumor.bam_index,
            normal_bam=realign_normal.bam_out,
            normal_bai=realign_normal.bam_index,
            reference_fa=reference_fa
    }

    output {
        File tumor_ig_minibam = realign_tumor.bam_out
        File tumor_ig_minibam_index = realign_tumor.bam_index
        File tumor_ig_minibam_idxstats = realign_tumor.bam_idxstats

        File normal_ig_minibam = realign_normal.bam_out
        File normal_ig_minibam_index = realign_normal.bam_index
        File normal_ig_minibam_idxstats = realign_normal.bam_idxstats

        File igcaller_filtered_tsv=IgCaller.igcaller_filtered_tsv
        File igcaller_IGH=IgCaller.igcaller_IGH
        File igcaller_IGK=IgCaller.igcaller_IGK
        File igcaller_IGL=IgCaller.igcaller_IGL
        File igcaller_CSR=IgCaller.igcaller_CSR
        File igcaller_IG_rearrangements=IgCaller.igcaller_IG_rearrangements
    }

}
