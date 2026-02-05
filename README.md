# Bipolar-englacial-bioinformatics-pipeline
This repository documents the bioinformatic tools, versions, and command-line arguments used in the analysis for O'Connor et al. 2026. ISME Communications. Bipolar investigation of near-surface glacial ice reveals an active microbial ecosystem driven by photosynthesis and chemolithoautotrophy.

---

No raw sequencing data are hosted in this repository. This repository serves as
a transparent record of the analytical workflow and software configuration.

---

## Overview of the workflow

The bioinformatic pipeline consisted of the following major steps:

1. Read quality control and adapter trimming  
2. Contaminant identification and removal using negative controls  
3. Metagenome assembly and taxonomic classification  
4. Genome binning, dereplication, and replication rate estimation  
5. Metatranscriptome decontamination and read mapping  
6. Functional and taxonomic annotation  

---

## Read preprocessing

Low-quality bases, reads, and adapter sequences were removed from both metagenomic
and metatranscriptomic datasets using **Trimmomatic v0.33** with the following parameters:

- `LEADING:3`
- `TRAILING:3`
- `SLIDINGWINDOW:4:15`

---

## Contaminant identification and removal

To identify contaminant sequences, metagenome reads from an artificial ice core
and a negative extraction control were co-assembled using **MEGAHIT v1.2.9**
with the `--meta-sensitive` setting.

This co-assembly was used to construct a contaminant reference database with
**DeconSeq v0.4.3**. DeconSeq was then used to remove contaminating sequences from
both glacier metagenomes and metatranscriptomes by filtering reads that mapped to
the negative control co-assembly.

---

## Metagenome assembly and taxonomic classification

Decontaminated metagenome reads were taxonomically classified using **Kaiju v1.9**
via the **KBase** platform.

Metagenome assembly was performed using **MEGAHIT v1.2.9**.

---

## Genome binning and refinement

Metagenome-assembled genomes (MAGs) were recovered using a combination of:

- **MetaBAT2 v2.15**
- **MaxBin2 v2.2.7**
- **SemiBin2 v1.5.1**

Resulting bins were dereplicated using **dRep v3.4.3**.
In situ replication rates of microbial populations represented by each MAG were
estimated using **iRep v1.10**.

---

## Metatranscriptome processing and read mapping

Ribosomal RNA sequences were removed from metatranscriptomic reads using
**SortMeRNA v4.3.6**.

Potential human DNA contamination was removed using the `removehuman` tool
from the **BBMap v38.92** package.

Filtered metatranscriptomic reads were aligned to the assembled metagenomes using
**Bowtie2 v2.5.1**, and transcript counts were generated using **HTSeq v2.0.2**.

---

## Functional and taxonomic annotation

Functional annotation of assembled metagenomes was performed by uploading assemblies
to the **JGI IMG/M** annotation pipeline.

---

## Software summary

| Tool | Version | Purpose |
|-----|--------|--------|
| Trimmomatic | 0.33 | Read trimming and QC |
| MEGAHIT | 1.2.9 | Metagenome assembly |
| DeconSeq | 0.4.3 | Contaminant removal |
| Kaiju | 1.9 | Taxonomic classification |
| MetaBAT2 | 2.15 | Genome binning |
| MaxBin2 | 2.2.7 | Genome binning |
| SemiBin2 | 1.5.1 | Genome binning |
| dRep | 3.4.3 | MAG dereplication |
| iRep | 1.10 | Replication rate estimation |
| SortMeRNA | 4.3.6 | rRNA removal |
| BBMap (removehuman) | 38.92 | Host contamination removal |
| Bowtie2 | 2.5.1 | Read alignment |
| HTSeq | 2.0.2 | Read counting |

---

## Notes

- Exact command-line arguments are documented in the `/scripts` directory.
- This repository is intended for transparency and reproducibility rather than
  full pipeline automation.
