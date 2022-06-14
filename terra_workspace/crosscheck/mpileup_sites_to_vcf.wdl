task mpileup {

    File bam
    File bam_index
    String prefix
    String prefix2
    String suffix = "vcf"
    File reference_fa
    File reference_fai
    File? sites_file
    String mpileup_args = ""
    String bctools_view_args = "-O v"

    String docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    Int? disk_size = ceil(size(bam, "GB") + size(reference_fa, "GB") + size(reference_fai, "GB"))
    Int? disk_size_buffer = "20"
    Int? disk_size_total = disk_size + disk_size_buffer
    Int? threads="1"
    Float? memory="3.75"
    Int? preempt = "4"

    command {
        set -euo pipefail

          mpileup_args_str=${mpileup_args}

          #Optionally add sites file
          if [ -n "${sites_file}" ]; then
              mpileup_args_str="$mpileup_args_str --targets-file ${sites_file}"
          fi

          #Run mpileup with flexible args and output format
          bcftools mpileup -f ${reference_fa} $mpileup_args_str ${bam} | bcftools view ${bctools_view_args} > "${prefix}${prefix2}.${suffix}"

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


workflow mpileup_workflow {
  call mpileup
}
