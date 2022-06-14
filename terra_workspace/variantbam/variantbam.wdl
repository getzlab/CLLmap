# BSD 3-Clause License
# Copyright (c) 2021, Broad Institute
# All rights reserved.
# Author: Binyamin A. Knisbacher

#Designed to work for both bam and cram input (output will be BAM)
task minibam {

    String sample_id
    File bam
    File? rule_file
    String? regions_opt = "-l" #-l gets mate too (see specs for other options: -g, -G, -L, -k)
    File? regions
	String? regions_str_opt = "-l"
    String? regions_str = ""
    String? command_line = ""
    String? minibam_descriptor = "mini"
    Int? region_pad = "0"

    String? docker = "gcr.io/broad-cga-bknisbac-wupo1/variantbam:1.0"
    Int? memory = "15"
    Int? disk_space_add = "50"
    Int? disk_space = ceil(size(bam ,"GB")) + disk_space_add
    Int? num_threads = "4"
    Int? num_preempt = "5"

    Int? num_threads_sort = num_threads - 1
    Int? memory_per_thread_sort = memory / num_threads

    command {
        set -euo pipefail

        #Define output prefix
        minibam_prefix=${sample_id}.${minibam_descriptor}

        # placeholders for optional outputs
        touch $minibam_prefix + ".bam.bai"
        touch $minibam_prefix + ".bam.idxstats"

        ### VariantBam command and parametres ###
        variant_args=" -v -b -t ${num_threads}"
        #Rule file for selection
        if [ -n "${rule_file}" ]; then
            variant_args="$variant_args --rules ${rule_file}"
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
        if [ -n "${region_pad}" && ${region_pad} > 0 ]; then
            variant_args="$variant_args -P ${region_pad}"
        fi
        #Free-style command line addition
        if [ -n "${command_line}" ]; then
            variant_args="$variant_args ${command_line}"
        fi
        #Run command
        variant ${bam} $variant_args > $minibam_prefix.bam


        ### sort + index (optional) ###
        samtools sort -O BAM -@ "${num_threads_sort}" -m ${memory_per_thread_sort}"G" -o $minibam_prefix.bam $minibam_prefix.bam
        samtools index $minibam_prefix.bam

        ### idxstats ###
        samtools idxstats $minibam_prefix.bam > $minibam_prefix.bam.idxstats
    }

    output {
        File minibam = "${sample_id}.${minibam_descriptor}.bam"
        File minibam_index = "${sample_id}.${minibam_descriptor}.bam.bai"
        File minibam_idxstats = "${sample_id}.${minibam_descriptor}.bam.idxstats"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Binyamin A. Knisbacher"
        email: "bknisbac@broadinstitute.org"
    }
}


workflow minibam_workflow {
    call minibam
}
