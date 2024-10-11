#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define default parameters
params.sra_ids = params.sra_ids ?: 'sra_ids.txt'  // File with SRA IDs, one per line
params.database_dir = params.database_dir ?: 'databases'
params.outdir = params.outdir ?: 'output'
params.keep_original_samples = params.keep_original_samples ?: false
params.help = params.help ?: false


// Show help message if requested
if (params.help) {
    log.info """
    Usage:
      nextflow run main.nf --sra_ids 'path/to/sra_ids.txt' --database_dir 'path/to/databases/' --outdir 'path/to/output/'

    Options:
      --sra_ids        : Path to file containing a list of SRA IDs (one per line).
      --database_dir   : Directory containing the required databases.
      --outdir         : Directory where output files will be saved.
      --help           : Display this help message.

    """
    exit 0
}

// Log the paths for confirmation
log.info """
    MetaChlam Pipeline
    ================================
    SRA IDs         : ${params.sra_ids}
    Database dir    : ${params.database_dir}
    Output dir      : ${params.outdir}
    """

// Process 1: Download data using fasterq-dump
process download_sample {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/sratools_conda_env.yml'

    tag "Download sample $sample_id"

    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")

    script:
    """
    # Download the sample using prefetch
    prefetch ${sample_id}

    # Create output directory
    mkdir -p ${params.outdir}/${sample_id}

    # Use fasterq-dump to split files
    fasterq-dump ${sample_id} --split-files -O ${params.outdir}/${sample_id}
    """
}


// strainscan         
process run_strainscan {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/strainscan_conda_env.yml'

    tag "filter $sample_id" 

    input: 
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("strainscan_output"), emit: strainscan_out

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        strainscan -i ${reads[0]} -j ${reads[1]} -d ${params.database_dir}/strainscan_db -o strainscan_output
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
    tuple val(sample_id), path("${sample_id}.classified*.fastq"), emit: classified
    tuple val(sample_id), path("${sample_id}.log"), emit: kraken_log

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        kraken2 --db ${params.database_dir}/LINtax_db --paired ${reads[0]} ${reads[1]} \
        --minimum-hit-groups 4 \
        --classified-out ${sample_id}.classified#.fastq >> ${sample_id}.log 2>&1

        kraken2 --db ${params.database_dir}/LINtax_db --paired ${reads[0]} ${reads[1]} \
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
    tuple val(sample_id), path("${sample_id}.LINreport.txt"), emit: lintax_out

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        python $baseDir/bin/report-lin.py --lin_file $baseDir/bin/lingroups.txt \
        --data_file ${params.database_dir}/LINtax_db/taxonomy/data.txt \
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
    tuple val(sample_id), path("${sample_id}.stats.tsv")
    tuple val(sample_id), path("${sample_id}.strains.tsv"), emit: strainge_out

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        straingst run -O -o ${sample_id} ${params.database_dir}/pan-genome-db_99.hdf5 ${sample_hdf5}
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
    tuple val(sample_id), path("${sample_id}.gather"), emit: sourmash_out

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        sourmash gather ${sample_sig} ${params.database_dir}/sourmash_db.sbt.zip -o ${sample_id}.gather
    """
}

// Ensure all processes synchronize outputs correctly
process run_global_report {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    
    input: 
    tuple val(sample_id),  path(kraken_log), path(strainscan_out), path(strainge_out), path(lintax_out), path(sourmash_out)

   // tuple val(sample_id), path(strainscan_out)
   // tuple val(sample_id), path(strainge_out)
   // tuple val(sample_id), path(lintax_out)
   // tuple val(sample_id), path(sourmash_out)
    
    output:
    path("${sample_id}_summary.csv"), emit: sample_summary

    script:
    """
        mkdir -p ${params.outdir}/${sample_id}
        python $baseDir/bin/global_report.py -kl ${kraken_log} -st ${strainge_out} \
        -lr ${lintax_out} -sc ${strainscan_out} \
        -sm ${sourmash_out} -s ${sample_id} \
        -o ${sample_id}_summary.csv
    """
}

// Final process: Cleanup original sample files (only retain output files)
process cleanup_sample_files {
    tag "Cleanup sample files $sample_id"

    input:
    tuple val(sample_id), path(original_files)

    when:
    !params.keep_original_samples // Only delete if the user hasn't requested to keep original files

    script:
    """
    rm -rf ${original_files[0]} ${original_files[1]}
    """
}

// Combine summaries process
process combine_summaries {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path collected_summaries

    output:
    path "Combined_summary.csv"

    script:
    """
    python -c "
import pandas as pd
summary_files = ['${collected_summaries.join("','")}']
combined_df = pd.concat([pd.read_csv(f) for f in summary_files], ignore_index=True)

# Ensure correct sorting by 'sample_id' (or replace 'sample_id' with your relevant column name)
if 'sample_id' in combined_df.columns:
    combined_df = combined_df.sort_values(by='sample_id')

combined_df.to_csv('Combined_summary.csv', index=False)
    "
    """
}

//workflow
workflow {
    // Define channel to read SRA IDs
    sra_ids_ch = Channel.fromPath(params.sra_ids).splitText().splitCsv()

    // Download samples using fasterq-dump
    download_ch = download_sample(sra_ids_ch)
    
    // Define channels for each process
    strainscan_out_ch = run_strainscan(download_ch)
    kraken_out_ch = run_kraken(download_ch)
    lintax_out_ch = run_lintax(kraken_out_ch.kraken_output, kraken_out_ch.kraken_report)
    strainge_kmer_out_ch = run_strainge_kmerize_sample(download_ch)
    strainge_out_ch = run_strainge(strainge_kmer_out_ch.sample_hdf5)
    sourmash_kmer_out_ch = run_sourmash_kmerize_sample(download_ch)
    sourmash_out_ch = run_sourmash(sourmash_kmer_out_ch.sample_sig)

   // Use `join` to join all the channels based on `sample_id`
    global_report_inputs = kraken_out_ch.kraken_log
                            .join(strainscan_out_ch.strainscan_out)
                            .join(strainge_out_ch.strainge_out)
                            .join(lintax_out_ch.lintax_out)
                            .join(sourmash_out_ch.sourmash_out)

    sample_summary_ch = run_global_report(global_report_inputs)

     // Cleanup original samples after the summary is generated
    cleanup_sample_files(download_ch)

    // Collect and combine all summary outputs
    sample_summary_ch
        .collect()
        .set { collected_summaries }

    combine_summaries(collected_summaries)
}