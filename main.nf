#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// set paths
params.reads = "$baseDir/sample/*_R{1,2}.fastq"
params.database_dir = "$baseDir/databases"
params.outdir = "$baseDir/output"
params.help = ""

def helpMessage() {
  log.info """
        Add Help Menu!!
        """ 
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// prints to the screen and to the log
log.info """
         MetaChlam (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

// strainscan         
process run_strainscan {

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("$sample_id")

    script:
    """
        strainscan -i ${reads[0]} -j ${reads[1]} -d $baseDir/databases/strainscan_db -o ${sample_id}/strainscan_output
    """
}

process run_LINtax {

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("$sample_id")

    script:
    """
        kraken2 -i ${reads[0]} -j ${reads[1]} -d $baseDir/databases/strainscan_db -o ${sample_id}/strainscan_output
    """
}

//workflow
workflow {
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    strainscan_out_ch = run_strainscan(reads)
}