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

1. Read Quality Control and Trimming
Tool: FastQC v.0.12.1
Forward and reverse reads were first checked with FastQC.

Tool: Trimmomatic v0.33
Applied to: Metagenomic and metatranscriptomic paired-end reads

java -jar trimmomatic-0.33.jar PE \
  -trimlog trimmomatic_log_file.txt \
  input_forward_read.fastq.gz input_reverse_read.fastq.gz \
  output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
  output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE-modified.fa:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36


Parameters:

Adapter removal: ILLUMINACLIP

Quality trimming: LEADING:3, TRAILING:3

Sliding window trimming: SLIDINGWINDOW:4:15

Minimum read length: MINLEN:36

2. Co-Assembly of Negative Controls

Tool: MEGAHIT v1.2.9

Negative control reads were co-assembled to generate a contaminant reference.

megahit \
  -1 output_forward_neg_sampl_1_paired.fq,output_forward_neg_sampl_2_paired.fq \
  -2 output_reverse_neg_sampl_1_paired.fq,output_reverse_neg_sampl_2_paired.fq \
  -o /output/folder

3. Contaminant Database Construction and Read Removal

Tool: DeconSeq

3.1 Create contaminant database
perl deconseq.pl \
  -f coassembled_negative_controls.fq \
  -dbs decontam_database

3.2 Remove contaminant reads
perl deconseq.pl \
  -f output_forward_paired.fq \
  -dbs decontam_database \
  -c 95 \
  -out_dir /output_directory

perl deconseq.pl \
  -f output_reverse_paired.fq \
  -dbs decontam_database \
  -c 95 \
  -out_dir /output_directory


Parameter:

-c 95 — minimum percent identity threshold

4. Metagenome Assembly

Tool: MEGAHIT v1.2.9

megahit \
  -1 decontaminated_forward_paired.fq \
  -2 decontaminated_reverse_paired.fq \
  -o /output/folder

5. Read Mapping to Assembly
5.1 Build Bowtie2 Index

Tool: Bowtie2 v2.5.1

bowtie2-build input_assembly.fa output_assembly_index_metagenome

5.2 Align Reads
bowtie2 \
  -x input_assembly_index_metagenome \
  -1 input_decontaminated_forward_reads.fq \
  -2 input_decontaminated_reverse_reads.fq \
  -S output_assembly_mapped.sam \
  --very-sensitive-local \
  --threads 60

5.3 Convert and Sort Alignments

Tool: samtools v1.16.1

samtools view -bS output_assembly_mapped.sam > output_assembly_mapped.bam
samtools sort output_assembly_mapped.bam -o output_assembly_mapped_sorted.bam


The sorted BAM file was used for genome binning.

6. Genome Binning

Three independent binning approaches were used.

6.1 MetaBAT2 (v2.15)
runMetaBat.sh \
  -m 2500 \
  assembly_contigs.fa \
  output_assembly_mapped_sorted.bam


Minimum contig length: 2500 bp.

6.2 MaxBin2 (v2.2.7)
perl /path/to/MaxBin-2.2.7/run_MaxBin.pl \
  -contig assembly_contigs.fa \
  -reads decontaminated_forward_paired.fq \
  -reads2 decontaminated_reverse_paired.fq \
  -out output_prefix \
  -thread 56

6.3 SemiBin2 (v1.5.1)
SemiBin2 single_easy_bin \
  --input-fasta assembly_contigs.fa \
  --input-bam output_assembly_mapped_sorted.bam \
  --environment global \
  --output output_prefix

7. MAG Quality Assessment

Tool: CheckM2 v1.0.0

checkm2 predict \
  -x fa \
  --input /path/to/input/bins \
  --output-directory /path/to/output/directory

8. Taxonomic Classification

Tool: GTDB-Tk v2.2.4

gtdbtk classify_wf \
  -x fa \
  --genome_dir /path/to/input/bins \
  --out_dir /path/to/output/directory

9. Dereplication

Tool: dRep v3.4.3

External completeness and contamination estimates generated with CheckM2
were supplied to dRep using the --genomeInfo flag.

dRep dereplicate /path/to/output/directory \
  -g /path/to/bins/*.fa \
  --genomeInfo checkM2_results.csv \
  --processors 56

10. In Situ Replication Rate Estimation

Tool: iRep v1.10

Bins used for iRep met the following criteria:

75% completeness

<2% contamination

<175 scaffolds per Mbp

10.1 Build Index
bowtie2-build input_bin.fa output_index_bin

10.2 Map Reads
bowtie2 \
  -x output_index_bin \
  -1 input_decontaminated_forward_reads.fq \
  -2 input_decontaminated_reverse_reads.fq \
  -S output_bin_mapped.sam \
  --very-sensitive-local \
  --threads 60

10.3 Run iRep
iRep \
  -f input_bin.fa \
  -s output_bin_mapped.sam \
  -o iRep_bin.iRep

11. Metatranscriptome Processing
11.1 Remove rRNA (SortMeRNA v4.3.6; DB v4.3)
sortmerna \
  --ref /path/to/database/smr_v4.3_default_db.fasta \
  --ref /path/to/database/smr_v4.3_fast_db.fasta \
  --ref /path/to/database/smr_v4.3_sensitive_db.fasta \
  --ref /path/to/database/smr_v4.3_sensitive_db_rfam_seeds.fasta \
  --reads /path/to/forward/trimmed/RNA/reads.fq.gz \
  --reads /path/to/reverse/trimmed/RNA/reads.fq.gz

11.2 Remove Human Sequences

BBMap v38.92 (removehuman)

./bbmap.sh \
  minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
  path=/path/to/removehuman-database/ qtrim=rl trimq=10 untrim -Xmx23g \
  in1=/path/to/mRNA/forward/reads.fq.gz \
  in2=/path/to/mRNA/reverse/reads.fq.gz \
  outu1=/path/to/forward/output/reads.fq.gz \
  outu2=/path/to/reverse/output/reads.fq.gz

12. Metatranscriptome Alignment and Gene Counting
12.1 Align Clean mRNA Reads
bowtie2-build input_contigs.fa output_index_contigs

bowtie2 \
  -x output_index_contigs \
  -1 /path/to/clean_forward_mRNA_reads.fq \
  -2 /path/to/clean_reverse_mRNA_reads.fq \
  -S output_mRNA_mapped.sam \
  --very-sensitive-local \
  --threads 60

12.2 GFF Identifier Reconciliation

Because JGI annotation replaces original contig identifiers, contig names were reconciled in R using the provided mapping file:

library(dplyr)

joined4 <- left_join(IC_MG_3_ga0599128_functional_annotation,
                     IC_MG_3_Ga0599128_contig_names_mapping,
                     by=c("V1" = "V2"))

joined4$V1 <- joined4$V1.y
joined5 <- joined4[,c(1:9)]

write.table(joined5,
            file="IC_MG_3_functional_annotation.gff",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)


Remove stray apostrophes:

sed -i "s/''//g" IC_MG_3_functional_annotation.gff

12.3 Sort Alignment File
samtools sort -o output_mRNA_mapped-sorted.sam output_mRNA_mapped.sam

12.4 Count Reads per Gene (HTSeq v2.0.2)
python -m HTSeq.scripts.count \
  -s no \
  -t CDS \
  -i ID \
  --nonunique=all \
  -r pos \
  -a 0 \
  /path/to/output_mRNA_mapped-sorted.sam \
  /path/to/IC_MG_3_functional_annotation.gff \
  > output_mRNA_counts.txt

Software Summary
Tool	Version	Purpose
Trimmomatic	0.33	Read trimming
MEGAHIT	1.2.9	Assembly
DeconSeq	—	Contaminant removal
Bowtie2	2.5.1	Alignment
samtools	1.16.1	BAM processing
MetaBAT2	2.15	Binning
MaxBin2	2.2.7	Binning
SemiBin2	1.5.1	Binning
CheckM2	1.0.0	MAG QC
GTDB-Tk	2.2.4	Taxonomy
dRep	3.4.3	Dereplication
iRep	1.10	Replication rate
SortMeRNA	4.3.6	rRNA removal
BBMap	38.92	Human sequence removal
HTSeq	2.0.2	Gene counting
Reproducibility Statement

This repository documents the analytical workflow, software versions, and parameters used in the study to support transparency and computational reproducibility. It is not intended to serve as a fully automated pipeline.
