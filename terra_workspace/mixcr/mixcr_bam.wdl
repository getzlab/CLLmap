workflow mixcr {
	call mixcr_task
}

task mixcr_task {
  File bam
  String id
  String sequence_type
  String? organism='hsa'
  String? receptor_type='xcr'
  String? other_args="--only-productive"
  String? docker="gcr.io/broad-cga-sanand-gtex/mixcr:latest"
  Int memoryGb
  Int diskSpaceGb
  Int num_preempt
  Int num_cores

	command {
        #placeholder for optional output
        touch ${id}.clns
        touch ${id}.contigs.clns

      java -Xmx${memoryGb}g -jar /opt/picard-tools/picard.jar SamToFastq \
          I=${bam} \
          FASTQ=fastq_R1.fastq \
          SECOND_END_FASTQ=fastq_R2.fastq

      java -Xmx${memoryGb}g -jar /opt/mixcr-3.0.10/mixcr.jar analyze shotgun \
          -s ${organism} \
          --starting-material ${sequence_type} \
          --align "-t ${num_cores}" \
          --assemble "-t ${num_cores}" \
          --report ${id}.report.txt \
          --receptor-type ${receptor_type} \
          ${other_args} \
          fastq_R1.fastq fastq_R2.fastq ${id}

        ls -l
        if [ ! -f "${id}.clonotypes.ALL.txt" ]; then
            cat ${id}.clonotypes.*.txt > ${id}.clonotypes.ALL.txt
        fi
        ls -l
	}

	output {
		File report = "${id}.report.txt"
		File mixcr_aligned = "${id}.vdjca"
		File mixcr_aligned_rescued_1 = "${id}.rescued_0.vdjca"
		File mixcr_aligned_rescued_2 = "${id}.rescued_1.vdjca"
		File mixcr_extended = "${id}.extended.vdjca"
		File mixcr_assembled = "${id}.clna"
        File mixcr_assembled_contigs = "${id}.clns"
        File mixcr_assembled_contigs2 = "${id}.contigs.clns"
    	File all_clones = "${id}.clonotypes.ALL.txt"
	}

	runtime {
		docker: "${docker}"
		memory: "${memoryGb} GB"
		cpu: "${num_cores}"
		disks: "local-disk ${diskSpaceGb} HDD"
   	preemptible: "${num_preempt}"
	}

  meta {
      author : "Shankara Anand" #modified by Binyamin Knisbacher
  }
}
