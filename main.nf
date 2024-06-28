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
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/strainscan_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("strainscan_output")

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        strainscan -i ${reads[0]} -j ${reads[1]} -d $baseDir/databases/strainscan_db -o strainscan_output
    """
}

process run_kraken {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/lintax_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.koutput"), emit: kraken_output
    tuple val(sample_id), path("${sample_id}.kreport"), emit: kraken_report

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        kraken2 --db $baseDir/databases/LINtax_db --paired ${reads[0]} ${reads[1]} \
        --minimum-hit-groups 4 \
        --confidence 0.45 \
        --output ${sample_id}.koutput \
        --report ${sample_id}.kreport
    """
}

process run_lintax {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/lintax_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(kraken_output)
    tuple val(sample_id), path(kraken_report)
    
    output:
    tuple val(sample_id), path("${sample_id}.LINreport.txt")

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        python $baseDir/bin/report-lin.py --lin_file $baseDir/bin/lingroups.txt \
        --data_file $baseDir/databases/LINtax_db/taxonomy/data.txt \
        --in_file_report ${kraken_report} \
        --in_file_output ${kraken_output} \
        --output ${sample_id}.LINreport.txt
    """
}

process run_strainge_kmerize_sample {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/strainge_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.hdf5"), emit: sample_hdf5

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        straingst kmerize -k 23 -o ${sample_id}.hdf5 ${reads[0]} ${reads[1]} 
    """
}

process run_strainge {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/strainge_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(sample_hdf5)
    
    output:
    tuple val(sample_id), path("${sample_id}.strainge.txt")

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        straingst run -o ${sample_id}.strainge.txt $baseDir/databases/pan-genome-db_99.hdf5 ${sample_hdf5}
    """
}

process run_sourmash_kmerize_sample {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/sourmash_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.zip"), emit: sample_sig

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        sourmash sketch dna -p scaled=1000,k=31 ${reads[0]} ${reads[1]} -o ${sample_id}.zip --name ${sample_id} 
    """
}

process run_sourmash {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/sourmash_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(sample_sig)
    
    output:
    tuple val(sample_id), path("${sample_id}.gather")

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        sourmash gather ${sample_sig} $baseDir/databases/sourmash_db.zip > ${sample_id}.gather
    """
}

//workflow
workflow {
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    strainscan_out_ch = run_strainscan(reads)
    kraken_out_ch = run_kraken(reads)
    lintax_out_ch = run_lintax(kraken_out_ch.kraken_output, kraken_out_ch.kraken_report)
    strainge_kmer_out_ch = run_strainge_kmerize_sample(reads)
    strainge_out_ch = run_strainge(strainge_kmer_out_ch.sample_hdf5)
    sourmash_kmer_out_ch = run_sourmash_kmerize_sample(reads)
    sourmash_out_ch = run_sourmash(sourmash_kmer_out_ch.sample_sig)
}
