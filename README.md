# MetaChlam
## Metagenomics tool for strain level identification of Chlamydia trachomatis strains
The MetaChlam pipeline is designed to classify and analyze metagenomic sequencing reads, specifically targeting Chlamydia trachomatis strains. The pipeline leverages multiple tools to process input reads, including StrainScan, LINtax, StrainGE, and Sourmash, and produces a combined summary report with results from each classifier.

## Requirements
Nextflow: The workflow manager.
Conda: Used to manage tool environments.

## Tools & Dependencies (handled installed by the pipeline)
This pipeline relies on several bioinformatics tools and packages:
- StrainScan
- Kraken2
- LINtax
- StrainGE
- Sourmash
- Python3

## Conda Environments (will be installed during the first run)
- strainscan_conda_env.yml: Environment for StrainScan.
- lintax_conda_env.yml: Environment for Kraken2 and LINtax.
- strainge_conda_env.yml: Environment for StrainGST.
- sourmash_conda_env.yml: Environment for Sourmash.

1. Install the pipeline
   ```
   git clone https://github.com/parul-sharma/MetaChlam.git
   ```

2. Set up the environment
   ```
   conda env create -f metachlam_env.yml
   ```

3. For executing the pipeline run:
   ```
   nextflow run main.nf
   ```

### Directory structure for running the MetaChLam pipeline:
MetaChlam/
├── bin/
│   ├── global_report.py           # Script to generate the global report
│   ├── lingroups.txt              # LIN group taxonomy file
│   └── report-lin.py              # Script to generate LIN reports
├── databases/
│   ├── LINtax_db/                 # LINtax database
│   ├── pan-genome-db_99.hdf5      # Pan-genome database for StrainGE
│   ├── sourmash_db.sbt.zip        # Sourmash database
│   └── strainscan_db/             # StrainScan database
├── envs/
│   ├── strainscan_conda_env.yml   # Conda environment for StrainScan
│   ├── lintax_conda_env.yml       # Conda environment for LINtax & Kraken2
│   ├── strainge_conda_env.yml     # Conda environment for StrainGST
│   └── sourmash_conda_env.yml     # Conda environment for Sourmash
├── main.nf                        # Main Nextflow pipeline
├── nextflow.config                # Nextflow configuration file
├── output/                        # Directory where all outputs are saved
│   └── (Generated output files)
└── sample/                        # Directory where all inputs are saved
    ├── sample_R1.fastq            # Example input read (R1)
    └── sample_R2.fastq            # Example input read (R2)


