workflow CrossCheckVCFvsBAM_Workflow {
    # WORKFLOW INPUT PARAMS
    # Pair Input
    # Sample tumor BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File tumorVCF
    # Sample normal BAM file (see https://samtools.github.io/hts-specs/SAMv1.pdf)
    File normalBam
    # Sample normal BAI file (an index for the normal BAM file) (see samtools index command http://www.htslib.org/doc/samtools.html)
    File normalBamIdx
    # Pair name, prefix for outputs
    String pairName
    # A list of common SNP loci to check for "fingerprint" mutations
    File HaplotypeDBForCrossCheck
    # Official GATK4 jar file to run the fingerprinting process with
    File GATK4_JAR
    # Preferred stringency level for the fingerprinting process
    String? validationStringencyLevel

    # RUNTIME OPTIONAL PARAMS
    String? diskGB_buffer
    String? preemptible
    String? memoryGB
    String? cpu

    # COMPUTE FILE SIZE
    Int gatk4_jar_size  = ceil(size(GATK4_JAR,  "G"))
    Int tumorVCF_size   = ceil(size(tumorVCF,   "G"))
    Int normalBam_size  = ceil(size(normalBam,  "G") + size(normalBamIdx,   "G"))

##############################

    # Program to check that all read groups within the set of BAM files appear to come from the same individual.
    call CrossCheckVCFvsBAM_Task {
        input:
            tumorVCF=tumorVCF,
            normalBam=normalBam,
            normalBamIdx=normalBamIdx,
            pairName=pairName,
            GATK4_JAR=GATK4_JAR,
            HaplotypeDBForCrossCheck=HaplotypeDBForCrossCheck,
            gatk4_jar_size=gatk4_jar_size,
            tumorVCF_size=tumorVCF_size,
            normalBam_size=normalBam_size,
            validationStringencyLevel=select_first([validationStringencyLevel, ""]),
            diskGB_buffer=select_first([diskGB_buffer, ""]),
            preemptible=select_first([preemptible, ""]),
            memoryGB=select_first([memoryGB, ""]),
            cpu=select_first([cpu, ""])
    }

    output {
        # Cross Check Lane Fingerprints Task
        File cross_check_fingprt_metrics=CrossCheckVCFvsBAM_Task.crossCheckMetrics
        File cross_check_fingprt_report=CrossCheckVCFvsBAM_Task.crossCheckReport
        Float cross_check_fingprt_min_lod_value=CrossCheckVCFvsBAM_Task.crossCheckMinLODValue
        String cross_check_fingprt_min_lod_lanes=CrossCheckVCFvsBAM_Task.crossCheckMinLODLanes
    }
}


task CrossCheckVCFvsBAM_Task {

    # TASK INPUT PARAMS
    File tumorVCF
    File normalBam
    File normalBamIdx
    String pairName
    File HaplotypeDBForCrossCheck
    String validationStringencyLevel
    File GATK4_JAR

    # FILE SIZE
    Int tumorVCF_size
    Int normalBam_size
    Int gatk4_jar_size

    # RUNTIME INPUT PARAMS
    String preemptible
    String diskGB_buffer
    String memoryGB
    String cpu

    # DEFAULT VALUES
    String default_cpu = "1"
    String default_memoryGB = "7"
    String default_preemptible = "1"
    String default_diskGB_buffer = "20"
    String default_stringencyLevel = "LENIENT"

    # COMPUTE MEMORY SIZE
    Int machine_memoryGB = if memoryGB != "" then memoryGB else default_memoryGB
    Int command_memoryGB = machine_memoryGB - 1
    Int machine_diskGB_buffer = if diskGB_buffer != "" then diskGB_buffer else default_diskGB_buffer

    # COMPUTE DISK SIZE
    Int diskGB = ceil(tumorVCF_size + normalBam_size + gatk4_jar_size + size(HaplotypeDBForCrossCheck, "G")
                    + machine_diskGB_buffer)

    String stringencyLevel = if validationStringencyLevel != "" then validationStringencyLevel else default_stringencyLevel

    parameter_meta {
        tumorVCF : "Sample tumor BAM file"
        normalBam : "Sample normal BAM file"
        normalBamIdx : "Sample normal BAI file (an index for the normal BAM file)"
        pairName : "Pair name, prefix for outputs"
        HaplotypeDBForCrossCheck : "A list of common SNP loci to check for \"fingerprint\" mutations"
        GATK4_JAR : "Official GATK4 jar file to run the fingerprinting process with"
        validationStringencyLevel : "Preferred stringency level for the fingerprinting process"
    }

    command <<<

        # e - exit when a command fails
        # u - exit when script tries to use undeclared variables
        # x - trace what gets executed
        # o - exit status of the last command that threw a non-zero exit code is returned
        set -euxo pipefail

        #Drop the seq entries which aren't in the BAM from haplotypeDB.
        PREPPED_HAPLOTYPE_DB=PreppedHaplotypeDB.txt
        echo -n "Number of snps in Haplotype DB before processing: " ; grep -vc '^@\|^#' ${HaplotypeDBForCrossCheck}
        /usr/local/bin/filter_not_in_bam_dict.pl ${normalBam} ${HaplotypeDBForCrossCheck} $PREPPED_HAPLOTYPE_DB
        echo -n "Number of snps in Haplotype DB after processing: " ; grep -vc '^@\|^#' $PREPPED_HAPLOTYPE_DB

        #CrosscheckLaneFingerprints[version=9]
        mkdir -v tmp
        java "-Xmx${command_memoryGB}g" -jar ${GATK4_JAR} CrosscheckFingerprints \
        -I ${tumorVCF} \
        -I ${normalBam} \
        -H $PREPPED_HAPLOTYPE_DB \
        --TMP_DIR `pwd`/tmp \
        --QUIET false \
        --EXIT_CODE_WHEN_MISMATCH 0 \
        --OUTPUT crosscheck.metrics \
        --VALIDATION_STRINGENCY ${stringencyLevel}

        #Produce crosscheck.stats.txt file for making the html report
        grep -v "#" crosscheck.metrics | sed 1d > crosscheck.metrics.stripped

        python /usr/local/bin/crosscheck_report.py crosscheck.metrics.stripped

        mv crosscheck.metrics ${pairName}.crosscheck.metrics
        mv report.html ${pairName}.report.html

    >>>

    runtime {
        docker: "gcr.io/broad-getzlab-workflows/cross_check_lane_fingerprints:1.0.0"
        preemptible: if preemptible != "" then preemptible else default_preemptible
        cpu: if cpu != "" then cpu else default_cpu
        disks: "local-disk ${diskGB} HDD"
        memory: machine_memoryGB + "GB"
    }

    output {
        File crossCheckMetrics="${pairName}.crosscheck.metrics"
        File crossCheckReport="${pairName}.report.html"
        Float crossCheckMinLODValue=read_float("crosscheck_min_lod_value.txt")
        String crossCheckMinLODLanes=read_string("crosscheck_min_lod_lanes.txt")
    }
}
