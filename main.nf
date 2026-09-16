#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads        = params.reads        ?: "$baseDir/sample/*_R{1,2}.fastq.gz"
params.database_dir = params.database_dir ?: "$baseDir/databases"
params.outdir       = params.outdir       ?: "$baseDir/output"
params.cleanup      = (params.cleanup != null) ? params.cleanup : true
params.help         = params.help ?: false

if (params.help) {
  log.info """
  Usage:
    nextflow run main.nf \\
      --reads 'path/*_R{1,2}.fastq.gz' \\
      --database_dir '/path/to/databases' \\
      --outdir '/path/to/output' \\
      [--cleanup true|false]

  Notes:
    - Kraken DB is Ct-only => classified reads are Ct reads
    - Ct_reads reported as paired reads: floor(classified/2)
  """
  exit 0
}

log.info """
MetaChlam (robust, mixed-infection aware)
===================================
reads        : ${params.reads}
database_dir : ${params.database_dir}
outdir       : ${params.outdir}
cleanup      : ${params.cleanup}
""".stripIndent()


// --------------------
// StrainScan
// --------------------
process run_strainscan {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/strainscan_conda_env.yml'
  tag "strainscan ${sample_id}"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("strainscan_output"), emit: strainscan_out

  script:
  """
  set +e
  mkdir -p strainscan_output

  strainscan -i ${reads[0]} -j ${reads[1]} -d ${params.database_dir}/strainscan_db -o strainscan_output

  status=\$?

  if [ \$status -ne 0 ]; then
      echo -e "sample_id\ttool\tstatus\tnote" > strainscan_output/strainscan_status.tsv
      echo -e "${sample_id}\tStrainScan\tfailed\tStrainScan crashed or produced no stable result; likely low-abundance/low-depth edge case" >> strainscan_output/strainscan_status.tsv
  else
      echo -e "sample_id\ttool\tstatus\tnote" > strainscan_output/strainscan_status.tsv
      echo -e "${sample_id}\tStrainScan\tsuccess\tcompleted" >> strainscan_output/strainscan_status.tsv
  fi

  exit 0
  """
}


// --------------------
// Kraken2 (Ct-only DB => classified reads are Ct reads; use log for classified counts)
// --------------------
process run_kraken {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/lintax_conda_env.yml'
  tag "kraken2 ${sample_id}"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.koutput"), emit: kraken_output
  tuple val(sample_id), path("${sample_id}.kreport"), emit: kraken_report
  tuple val(sample_id), path("${sample_id}.classified*.fastq.gz"), emit: classified
  tuple val(sample_id), path("${sample_id}.extract_read.log"), emit: extract_read_log

  script:
  """
  set -euo pipefail
  
  # Detect whether inputs are gzipped
  KRAKEN_GZIP=""
  if [[ "${reads[0]}" == *.gz ]]; then
    KRAKEN_GZIP="--gzip-compressed"
  fi

  # Run 1: produce classified reads + log (Ct-only DB => classified == Ct reads)
  kraken2 --db ${params.database_dir}/kraken_db/ct_db --paired \$KRAKEN_GZIP ${reads[0]} ${reads[1]} \\
    --minimum-hit-groups 4 \\
    --confidence 0.3 \\
    --output ${sample_id}.ctscreen.tmp.koutput > /dev/null 2> ${sample_id}.ctscreen.tmp.stderr

  ${baseDir}/bin/extract_kraken_reads.py \\
    -k ${sample_id}.ctscreen.tmp.koutput \\
    -s1 ${reads[0]} \\
    -s2 ${reads[1]} \\
    -t 813 \\
    --fastq-output \\
    -o ${sample_id}.classified_1.fastq \\
    -o2 ${sample_id}.classified_2.fastq > ${sample_id}.extract_read.log 2>&1

  rm -f ${sample_id}.ctscreen.tmp.koutput ${sample_id}.ctscreen.tmp.stderr


  # Compress classified fastqs (keeps disk usage down)
  if ls ${sample_id}.classified*.fastq >/dev/null 2>&1; then
    gzip -f ${sample_id}.classified*.fastq
  fi

  # Run 2: produce koutput/kreport
  kraken2 --db ${params.database_dir}/LINtax_db --paired \$KRAKEN_GZIP ${reads[0]} ${reads[1]} \\
    --minimum-hit-groups 4 \\
    --confidence 0.45 \\
    --output ${sample_id}.koutput \\
    --report ${sample_id}.kreport
  """
}



// --------------------
// LINtax (stringent)
/// --------------------
process run_lintax {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/lintax_conda_env.yml'
  tag "lintax ${sample_id}"

  input:
  tuple val(sample_id), path(kraken_output)
  tuple val(sample_id), path(kraken_report)

  output:
  tuple val(sample_id), path("${sample_id}.LINreport.txt"), emit: lintax_out

  script:
  """
  set -euo pipefail
  python $baseDir/bin/report-lin.py --lin_file $baseDir/bin/lingroups.txt \\
    --data_file ${params.database_dir}/LINtax_db/taxonomy/data.txt \\
    --in_file_report ${kraken_report} \\
    --in_file_output ${kraken_output} \\
    --output ${sample_id}.LINreport.txt

  [ -f ${sample_id}.LINreport.txt ] || touch ${sample_id}.LINreport.txt
  """
}


// --------------------
// StrainGE kmerize (supports gz)
// --------------------
process run_strainge_kmerize_sample {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/strainge_conda_env.yml'
  tag "straingst-kmerize ${sample_id}"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.hdf5"), emit: sample_hdf5

  script:
  """
  set -euo pipefail
  straingst kmerize -k 23 -o ${sample_id}.hdf5 ${reads[0]} ${reads[1]}
  """
}


// --------------------
// StrainGE run (robust to 0 Ct reads crash)
// Writes a placeholder strains.tsv so reporting never misses sample
// --------------------
process run_strainge {
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda 'envs/strainge_conda_env.yml'
    tag "straingst-run ${sample_id}"

    input:
    tuple val(sample_id), path(sample_hdf5)

    output:
    tuple val(sample_id), path("${sample_id}.stats.tsv")
    tuple val(sample_id), path("${sample_id}.strains.tsv"), emit: strainge_out

    script:
    """
    set -u

    make_placeholders () {
      # Minimal valid placeholders so downstream parsing never breaks
      cat > ${sample_id}.strains.tsv <<'EOF'
strain\tgcov
NA\t0
EOF

      cat > ${sample_id}.stats.tsv <<'EOF'
metric\tvalue
status\tno_ct_reads_or_strainge_failed
EOF
    }

    # If the sample hdf5 is missing/empty, just write placeholders
    if [ ! -s "${sample_hdf5}" ]; then
      echo "INFO: ${sample_id} sample_hdf5 missing/empty -> placeholders" >&2
      make_placeholders
      exit 0
    fi

    # Run StrainGE; if it fails (e.g. NaN crash), do not fail the pipeline
    if ! straingst run -O -o ${sample_id} $baseDir/databases/pan-genome-db_99.hdf5 "${sample_hdf5}" ; then
      echo "INFO: ${sample_id} StrainGE failed -> placeholders" >&2
      make_placeholders
      exit 0
    fi

    # If outputs weren’t produced for any reason, still create placeholders
    if [ ! -s ${sample_id}.strains.tsv ] || [ ! -s ${sample_id}.stats.tsv ]; then
      echo "INFO: ${sample_id} StrainGE outputs missing/empty -> placeholders" >&2
      make_placeholders
    fi
    """
}




// --------------------
// Sourmash
// --------------------
process run_sourmash_kmerize_sample {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/sourmash_conda_env.yml'
  tag "sourmash-sketch ${sample_id}"

  input:
  tuple val(sample_id), path(reads)

  output:
  tuple val(sample_id), path("${sample_id}.zip"), emit: sample_sig

  script:
  """
  set -euo pipefail
  sourmash sketch dna -p scaled=1000,k=31 ${reads[0]} ${reads[1]} -o ${sample_id}.zip --name ${sample_id}
  """
}

process run_sourmash {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  conda 'envs/sourmash_conda_env.yml'
  tag "sourmash-gather ${sample_id}"

  input:
  tuple val(sample_id), path(sample_sig)

  output:
  tuple val(sample_id), path("${sample_id}.gather"), emit: sourmash_out

  script:
  """
  set -euo pipefail
  sourmash gather ${sample_sig} ${params.database_dir}/sourmash_db.sbt.zip -o ${sample_id}.gather
  [ -f ${sample_id}.gather ] || touch ${sample_id}.gather
  """
}


// --------------------
// Global report: writes BOTH summary and evidence
// --------------------
process run_global_report {
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'
  tag "global-report ${sample_id}"

  input:
  tuple val(sample_id),
        path(extract_read_log),
        path(strainscan_out),
        path(strainge_out),
        path(lintax_out),
        path(sourmash_out)

  output:
  tuple val(sample_id), path("${sample_id}_summary.csv"), emit: sample_summary
  tuple val(sample_id), path("${sample_id}_evidence.tsv"), emit: sample_evidence

  script:
  """
  set -euo pipefail
  python $baseDir/bin/global_report.py \\
    -el ${extract_read_log} \\
    -st ${strainge_out} \\
    -lr ${lintax_out} \\
    -sm ${sourmash_out} \\
    -sc ${strainscan_out} \\
    -s ${sample_id} \\
    -o ${sample_id}_summary.csv \\
    -e ${sample_id}_evidence.tsv
  """
}


// --------------------
// Cleanup big intermediates after report exists
// --------------------
process cleanup_big_files {
  tag "cleanup ${sample_id}"

  input:
  tuple val(sample_id),
        path(sample_summary),
        path(sample_evidence),
        path(koutput),
        path(extract_read_log),
        path(sample_hdf5)

  when:
  params.cleanup

  script:
  """
  set -euo pipefail
  rm -f "${params.outdir}/${sample_id}/${sample_id}.koutput" || true
  #rm -f "${params.outdir}/${sample_id}/${sample_id}.extract_read.log"    || true
  rm -f "${params.outdir}/${sample_id}/${sample_id}.hdf5"   || true

  # staged copies
  rm -f "${koutput}"  "${sample_hdf5}" || true
  """
}


// --------------------
// Combine summaries + evidence
// --------------------
process combine_reports {
  publishDir "${params.outdir}", mode: 'copy'
  tag "combine-reports"

  input:
  path collected_summaries
  path collected_evidence

  output:
  path "Combined_summary.csv"
  path "Combined_evidence.tsv"

  script:
  """
  python - <<'PY'
import pandas as pd

# summaries
#summ_files = ${collected_summaries.collect { "'" + it.toString() + "'" }}
#summ = pd.concat([pd.read_csv(f) for f in summ_files], ignore_index=True)
#if 'Sample' in summ.columns:
#    summ = summ.sort_values(by='Sample')
#summ.to_csv("Combined_summary.csv", index=False)

# evidence
#ev_files = ${collected_evidence.collect { "'" + it.toString() + "'" }}
#ev = pd.concat([pd.read_csv(f, sep="\\t") for f in ev_files], ignore_index=True)
#if 'Sample' in ev.columns:
#    ev = ev.sort_values(by=['Sample','Tool','Hit_rank'])
#ev.to_csv("Combined_evidence.tsv", sep="\\t", index=False)
#PY

# summaries
summ_files = ${collected_summaries.collect { "'" + it.toString() + "'" }}
summ_dfs = []
for f in summ_files:
    try:
        df = pd.read_csv(f)
        if not df.empty:
            summ_dfs.append(df)
    except Exception as e:
        print(f"[WARN] Could not read summary {f}: {e}")

summ = pd.concat(summ_dfs, ignore_index=True) if summ_dfs else pd.DataFrame()
if 'Sample' in summ.columns:
    summ = summ.sort_values(by='Sample')
summ.to_csv("Combined_summary.csv", index=False)

# evidence
ev_files = ${collected_evidence.collect { "'" + it.toString() + "'" }}
ev_dfs = []
for f in ev_files:
    try:
        df = pd.read_csv(f, sep="\\t")
        if not df.empty:
            ev_dfs.append(df)
    except Exception as e:
        print(f"[WARN] Could not read evidence {f}: {e}")

ev = pd.concat(ev_dfs, ignore_index=True) if ev_dfs else pd.DataFrame()
if {'Sample','Tool','Hit_rank'}.issubset(ev.columns):
    ev = ev.sort_values(by=['Sample','Tool','Hit_rank'])
ev.to_csv("Combined_evidence.tsv", sep="\\t", index=False)
PY
  """
}


// --------------------
// Workflow
// --------------------
workflow {
  reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

  strainscan_out_ch     = run_strainscan(reads_ch)
  kraken_out_ch         = run_kraken(reads_ch)
  lintax_out_ch         = run_lintax(kraken_out_ch.kraken_output, kraken_out_ch.kraken_report)
  strainge_kmer_out_ch  = run_strainge_kmerize_sample(reads_ch)
  strainge_out_ch       = run_strainge(strainge_kmer_out_ch.sample_hdf5)
  sourmash_kmer_out_ch  = run_sourmash_kmerize_sample(reads_ch)
  sourmash_out_ch       = run_sourmash(sourmash_kmer_out_ch.sample_sig)

  // Join inputs by sample_id; since StrainGE writes a placeholder, samples won't disappear here.
  global_inputs =
    kraken_out_ch.extract_read_log
      .join(strainscan_out_ch.strainscan_out)
      .join(strainge_out_ch.strainge_out)
      .join(lintax_out_ch.lintax_out)
      .join(sourmash_out_ch.sourmash_out)

  report_ch = run_global_report(global_inputs)

  // cleanup after report exists
  cleanup_inputs =
    report_ch.sample_summary
      .join(report_ch.sample_evidence)
      .join(kraken_out_ch.kraken_output)
      .join(kraken_out_ch.extract_read_log)
      .join(strainge_kmer_out_ch.sample_hdf5)

  cleanup_big_files(cleanup_inputs)

  // collect for combined outputs
  report_ch.sample_summary
    .map { sid, f -> f }
    .collect()
    .set { collected_summaries }

  report_ch.sample_evidence
    .map { sid, f -> f }
    .collect()
    .set { collected_evidence }

  combine_reports(collected_summaries, collected_evidence)
}
