task extract_fingerprint {

    File bam
    File bam_index
    String prefix
    String prefix2
    String suffix = "vcf"
    File reference_fa
    File reference_fai
    File reference_dict
    File haplotype_map
    String extract_fingerprint_args = ""

    String docker = "broadinstitute/gatk:4.2.0.0"
    Int? disk_size = ceil(size(bam, "GB") + size(bam_index, "GB") + size(reference_fa, "GB") + size(reference_fai, "GB") + size(haplotype_map, "GB"))
    Int? disk_size_buffer = "20"
    Int? disk_size_total = disk_size + disk_size_buffer
    Int? threads="1"
    Float? memory="3.75"
    Int? preempt = "4"

    command {
        set -euo pipefail

          #Run ExtractFingerprint with flexible args
          gatk ExtractFingerprint -I ${bam} -H ${haplotype_map} -R ${reference_fa} -O "${prefix}${prefix2}.${suffix}" ${extract_fingerprint_args}

    }
    output {
        File vcf="${prefix}${prefix2}.${suffix}"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_size_total} HDD"
        cpu: "${threads}"
        preemptible: "${preempt}"
    }

    meta {
        author: "Binyamin A. Knisbacher"
    }
}


workflow extract_fingerprint_workflow {
  call extract_fingerprint
}
