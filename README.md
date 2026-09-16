# MetaChlam

## Strain-level detection and typing of *Chlamydia trachomatis* from metagenomic sequencing data

📄 **Preprint:** *Detecting and typing Chlamydia trachomatis strains in metagenomes using the MetaChlam pipeline*  
Sharma P, Dean D, Read TD. bioRxiv (2026).  
DOI: https://doi.org/10.64898/2026.08.18.745514 


MetaChlam is a **Nextflow-based** bioinformatics pipeline for detecting and characterizing *Chlamydia trachomatis* directly from metagenomic sequencing data. Rather than relying on a single classifier, MetaChlam integrates multiple complementary strain typing approaches into a single automated workflow and produces a consensus report summarizing strain assignments across multiple methods.

Although this repository is currently distributed with **custom databases for *Chlamydia trachomatis***, the workflow itself is modular and can be adapted for other bacterial pathogens by replacing the organism-specific reference databases used by each classification tool.

---

## Features

- Nextflow-based reproducible workflow
- Automatic Conda environment management
- Supports paired-end Illumina metagenomic reads
- Quick validation of input FASTQ files before analysis
- Integrates four complementary strain typing methods:
  - LINtax
  - StrainScan
  - StrainGE
  - Sourmash
- Generates a unified summary report combining results from all tools
- Includes curated *C. trachomatis* reference databases
- Easily adaptable to other bacterial species by substituting custom databases

---

## Pipeline overview

```text
Paired-end FASTQ files
        │
        ▼
Input validation
        │
        ▼
┌──────────────┬──────────────┬──────────────┬──────────────┐
│   LINtax     │  StrainScan  │   StrainGE   │   Sourmash   │
└──────────────┴──────────────┴──────────────┴──────────────┘
        │
        ▼
Consensus integration
        │
        ▼
Global MetaChlam report
```

---

## Requirements

- Linux
- Nextflow (≥24 recommended)
- Conda or Mamba

The individual software packages are installed automatically during the first execution.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/parul-sharma/MetaChlam.git
cd MetaChlam
```

Create the MetaChlam environment:

```bash
conda env create -f metachlam_env.yml
conda activate metachlam
```

The remaining Conda environments required by individual tools will be created automatically when the pipeline is first run.

---

## Running MetaChlam

Example:

```bash
nextflow run main.nf \
    --reads "/path/to/*_{R1,R2}.fastq.gz" \
    --outdir results
```

Cluster-specific execution settings can be modified in `nextflow.config`.

---

## Pipeline components

| Tool | Purpose |
|------|---------|
| LINtax | Genome similarity-based lineage classification |
| StrainScan | Marker-based strain identification |
| StrainGE | Strain detection in low-abundance metagenomes |
| Sourmash | MinHash-based genome similarity search |

---

## Repository structure

```text
MetaChlam/
├── bin/                  Helper scripts
├── databases/            Organism-specific reference databases
├── envs/                 Conda environment files
├── main.nf               Main Nextflow workflow
├── nextflow.config       Pipeline configuration
├── conf/                 Optional execution profiles
└── README.md
```

---

## Output

MetaChlam generates:

- Individual reports from each strain typing method
- A combined consensus summary report
- Tool-specific log files
- Nextflow execution reports and logs

---

## Customizing MetaChlam

The workflow is currently configured with curated ***Chlamydia trachomatis*** databases for LINtax, StrainScan, StrainGE, and Sourmash.

To adapt MetaChlam for another bacterial species, replace the organism-specific reference databases for each tool with databases built for the target organism. No changes to the overall workflow architecture are required.

For guidance on adapting MetaChlam to additional pathogens, please contact:

**Parul Sharma**  
Department of Infectious Diseases  
Emory University  
📧 parul.sharma2@emory.edu

---

## Citation

If you use MetaChlam in your research, please cite:

> Sharma P, Dean D, Read TD. 
> *Detecting and typing Chlamydia trachomatis strains in metagenomes using the MetaChlam pipeline* bioRxiv (2026).   
DOI: https://doi.org/10.64898/2026.08.18.745514 

---

## License

MIT License