# **Africonus crotchii** genome assembly
 This repository contains the computational pipeline, scripts, and genomic resources for the haplotype-resolved, chromosome-level genome assembly of the Cabo Verde venomous snail *Africonus crotchii*.
 
 ## Project Overview
 Bioinformatic pipeline for haplotype-resolved genome assembly, annotation, and comparative genomic analysis of *Africonus crotchii* using PacBio HiFi, Omni-C, and RNA-seq data.
 The workflow generates and evaluates two haplotype-resolved chromosome-scale genome assemblies, each organized into 35 pseudochromosomes.

 ## Genome assembly highlights
 We combined PacBio HiFi (CCS) and Omni-C sequencing technologies to generate two haplotype-resolved, chromosome-scale genome assemblies of *Africonus crotchii*.
 | Assembly statistics | Haplotype 1 | Haplotype 2 |
 |---|---:|---:|
 | **Assembly size** | 2.28 Gb | 3.18 Gb |
 | **Scaffold N50** | 60.23 Mb | 78.58 Mb |
 | **BUSCO completeness** (`mollusca_odb12`) | 75.7% | 91.6% |
 | **Pseudochromosomes** | 35 | 35 |
 | **Genome anchored to pseudochromosomes** | 86.86% | 85.23% |

## Bioinformatic pipeline

The complete command-line workflow is available in
[the Bash pipeline](scripts)

### 1. Genome assembly

**Software:** Hifiasm, BLAST, BlobTools.

- **Input:** PacBio HiFi (CCS) reads.
- **Workflow:** *De novo* assembly using Hifiasm, followed by contamination screening against the NCBI nucleotide database and taxonomic evaluation with BlobTools.
- **Output:** Haplotype-resolved contig assemblies.

### 2. Chromosome-level scaffolding

**Software:** BWA, HiRise.

- **Input:** Haplotype-specific contig assemblies and Omni-C paired-end reads.
- **Workflow:** Independent scaffolding of both haplotypes using Omni-C proximity-ligation data. HiRise was used to identify and correct misjoins and generate chromosome-scale scaffolds.
- **Output:** Two chromosome-scale assemblies, each comprising 35 pseudochromosomes.

*Note: HiRise scaffolding was performed externally using proprietary software.*

### 3. Genome quality assessment

**Software:** QUAST, BUSCO, Merqury, Meryl, minimap2, SAMtools.

- **Input:** Final haplotype assemblies and PacBio HiFi reads.
- **Workflow:** Assessment of genome contiguity, completeness, consensus accuracy, and read-mapping statistics. BUSCO was run using the Metazoa and Mollusca datasets, while Merqury evaluated each haplotype independently and their combined k-mer completeness.
- **Output:** Assembly statistics, BUSCO completeness scores, consensus quality values (QV), and k-mer completeness estimates.

### 4. Repeat annotation

**Software:** RepeatModeler2, RepeatMasker.

- **Input:** Final haplotype assemblies.
- **Workflow:** Haplotype-specific *de novo* repeat identification, classification, and soft-masking. Repeat landscapes were generated using RepeatMasker utilities to investigate transposable element abundance and sequence divergence.
- **Output:** Soft-masked genome assemblies, repeat annotations, and repeat landscapes.

### 5. Gene prediction

**Software:** Trimmomatic, STAR, SAMtools, BRAKER3.

- **Input:** Soft-masked haplotype assemblies, venom gland RNA-seq reads, and protein datasets from related cone snails and Metazoa OrthoDB v12.
- **Workflow:** RNA-seq preprocessing, read alignment, and evidence-supported gene prediction combining transcriptomic and protein information.
- **Output:** Haplotype-specific genome annotations (GFF3) and predicted gene models.

### 6. Haplotype comparison

**Software:** minimap2, D-GENIES, MCScanX, SynVisio, BUSCO, and custom Python scripts.

- **Input:** Both haplotype assemblies, BRAKER3 annotations, repeat annotations, BUSCO results, and reference cone snail genomes.
- **Workflow:**
  - Comparison of gene, exon, and intron length distributions.
  - Comparison of repeat landscapes.
  - Whole-genome alignments and structural collinearity.
  - Ortholog-level BUSCO completeness comparisons and investigation of missing genomic regions.
  - Interspecific synteny analyses with *L. ventricosus* and *K. canariensis*.
- **Output:** Comparative statistics, boxplots, dot plots, BUSCO comparison tables, and chromosomal synteny visualizations.

## Software and dependencies

The pipeline combines Bash command-line tools, Perl utilities, and custom Python scripts. Software versions and relevant parameters are documented in the accompanying pipeline script.

Analyses were performed independently for each haplotype unless otherwise indicated. Steps involving proprietary software or online visualization platforms are documented but cannot be reproduced directly from the command-line script.

## Data availability

The two haplotype-resolved genome assemblies of *Africonus crotchii* have been deposited in DDBJ/ENA/GenBank.

### Genome Assemblies & Accessions

| | Haplotype 1 | Haplotype 2 |
|---|---|---|
| **WGS project** | JBSYHA000000000 | JBSYHB000000000 |
| **Assembly version 01** | JBSYHA010000000 | JBSYHB010000000 |
| **BioProject** | [PRJNA1381619](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1381619) | [PRJNA1381620](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1381620) |
| **BioSample** | SAMN53902219 | SAMN53902219 |

### Raw Sequencing Data (SRA Accessions)

The genome assemblies were generated using PacBio HiFi sequencing and Omni-C proximity-ligation data. Venom gland RNA-seq data were used to support gene prediction.

| **Platform** | **Data type** | **SRA accession** |
|---|---|---|
| PacBio Sequel II | HiFi (CCS) long reads | SRR36508567 |
| Illumina HiSeq X | Omni-C proximity-ligation reads | SRR36508566 |
| Illumina | Venom gland RNA-seq | SRR36464808 - SRR36464812 |

## Citation

If you use the data, genome assembly, or scripts from this repository, please cite:

Herráez-Pérez, A., Pardos-Blas, J.R., & Zardoya, R. 2026. Haplotype-Resolved Chromosome-Level Genome of the Cabo Verde Venomous Snail *Africonus crotchii*. *Submitted to G3: Genes, Genomes, Genetics

## License

See the accompanying LICENSE file for licensing terms.
