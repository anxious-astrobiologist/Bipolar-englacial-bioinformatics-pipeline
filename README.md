# Bipolar-englacial-bioinformatics-workflow
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

The majority of tools were installed in a Conda environment.

---

### 1. Quality Control

Tool: Trimmomatic v0.33
Performed on: Metagenomic and metatranscriptomic reads
```
java -jar trimmomatic-0.33.jar PE -trimlog trimmomatic_log_file.txt \
input_forward_read.fastq.gz input_reverse_read.fastq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE-modified.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
### 2. Decontamination
#### 2.1 Co-assembly of Negative Controls

Tool: MegaHit v1.2.9
```
megahit -1 output_forward_neg_sampl_1_paired.fq,output_forward_neg_sampl_2_paired.fq \
-2 output_reverse_neg_sampl_1_paired.fq,output_reverse_neg_sampl_2_paired.fq \
-o /output/folder
```
#### 2.2 Create Decontamination Database

Tool: DeconSeq v0.4.3
```
perl deconseq.pl -f coassembled_negative_controls.fq -dbs decontam_database
```
#### 2.3 Remove Contaminant Reads
```
perl deconseq.pl -f output_forward_paired.fq -dbs decontam_database -c 95 -out_dir /output_directory
perl deconseq.pl -f output_reverse_paired.fq -dbs decontam_database -c 95 -out_dir /output_directory
```

-c 95 = 95% identity cutoff for contaminant matches

Taxonomic classification of the decontaminated reads was performed using Kaiju v1.9 through the KBase platform.

### 3. Metagenome Assembly & Annotation
#### 3.1 Metagenome Assembly
Tool: MegaHit v1.2.9
```
megahit -1 decontaminated_forward_paired.fq \
-2 decontaminated_reverse_paired.fq \
-o /output/folder
```
#### 3.2 Metagenome Functional Annotation

Metagenomes were annotated using the JGI IMG/M pipeline.

### 4. Read Mapping to Assembly
#### 4.1 Build Bowtie2 Index

Tool: Bowtie2 v2.5.1
```
bowtie2-build input_assembly.fa output_assembly_index
```
#### 4.2 Align Reads
```
bowtie2 -x output_assembly_index \
-1 input_decontaminated_forward_reads.fq \
-2 input_decontaminated_reverse_reads.fq \
-S output_assembly_mapped.sam \
--very-sensitive-local --threads 60
```
#### 4.3 Convert and Sort

Tool: samtools v1.16.1
```
samtools view -bS output_assembly_mapped.sam > output_assembly_mapped.bam
samtools sort output_assembly_mapped.bam -o output_assembly_mapped_sorted.bam
```
### 5. Genome Binning
MetaBat2 v2.15
```
runMetaBat.sh -m 2500 assembly_contigs.fa output_assembly_mapped_sorted.bam
```
MaxBin2 v2.2.7
```
perl run_MaxBin.pl -contig assembly_contigs.fa \
-reads decontaminated_forward_paired.fq \
-reads2 decontaminated_reverse_paired.fq \
-out output_prefix -thread 56
```
SemiBin2 v1.5.1
```
SemiBin2 single_easy_bin --input-fasta assembly_contigs.fa \
--input-bam output_assembly_mapped_sorted.bam \
--environment global \
--output output_prefix
```
### 6. Bin Quality Assessment

Tool: CheckM2 v1.0.0
```
checkm2 predict -x fa \
--input /path/to/input/bins \
--output-directory /path/to/output/directory
```
### 7. Bin Taxonomic Classification

Tool: GTDB-Tk v2.2.4
```
gtdbtk classify_wf -x fa \
--genome_dir /path/to/input/bins \
--out_dir /path/to/output/directory
```
### 8. Dereplication

Tool: dRep v3.4.3

Used CheckM2-derived completeness/contamination values via --genomeInfo.
```
dRep dereplicate /path/to/output/directory \
-g /path/to/bins/*.fa \
--genomeInfo checkM2_results.csv \
--processors 56
```
### 9. In Situ Replication Rate

Tool: iRep v1.10

Build Index
```
bowtie2-build input_bin.fa output_index_bin
```
Map Reads
```
bowtie2 -x output_index_bin \
-1 input_decontaminated_forward_reads.fq \
-2 input_decontaminated_reverse_reads.fq \
-S output_bin_mapped.sam \
--very-sensitive-local --threads 60
```
Run iRep
```
iRep -f input_bin.fa -s output_bin_mapped.sam -o iRep_bin.iRep
```

Note: Bins must be >75% complete, <2% contamination, and <175 scaffolds/Mbp.

### 10. Metatranscriptome Processing
#### 10.1 Remove rRNA

Tool: SortMeRNA v4.3.6
```
sortmerna \
--ref /path/to/database/smr_v4.3_default_db.fasta \
--ref /path/to/database/smr_v4.3_fast_db.fasta \
--ref /path/to/database/smr_v4.3_sensitive_db.fasta \
--ref /path/to/database/smr_v4.3_sensitive_db_rfam_seeds.fasta \
--reads forward_reads.fq.gz \
--reads reverse_reads.fq.gz
```
#### 10.2 Remove Human Sequences

Tool: BBMap v38.92 (removehuman)
```
./bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
path=/path/to/removehuman-database/ qtrim=rl trimq=10 untrim -Xmx23g \
in1=forward_reads.fq.gz \
in2=reverse_reads.fq.gz \
outu1=clean_forward_reads.fq.gz \
outu2=clean_reverse_reads.fq.gz
```
### 11. Gene Expression Quantification
#### 11.1 Align mRNA to Metagenome
```
bowtie2-build input_contigs.fa output_index_contigs

bowtie2 -x output_index_contigs \
-1 clean_forward_mRNA_reads.fq \
-2 clean_reverse_mRNA_reads.fq \
-S output_mRNA_mapped.sam \
--very-sensitive-local --threads 60
```
#### 11.2 Prepare GFF File (R)

Contig names were corrected using JGI-provided mapping files in R.
```
library(dplyr)

joined4 <- left_join(annotation_file, contig_mapping, by=c("V1" = "V2"))
joined4$V1 <- joined4$V1.y
joined5 <- joined4[,c(1:9)]

write.table(joined5,
            file="functional_annotation.gff",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)
```

Remove apostrophes (from a terminal):
```
sed -i "s/''//g" functional_annotation.gff
```
#### 11.3 Sort Alignment File
```
samtools sort -o output_mRNA_mapped-sorted.sam output_mRNA_mapped.sam
```
#### 11.4 Count Reads Per Gene

Tool: HTSeq v2.0.2
```
python -m HTSeq.scripts.count \
-s no -t CDS -i ID --nonunique=all -r pos -a 0 \
output_mRNA_mapped-sorted.sam \
functional_annotation.gff \
> output_mRNA_counts.txt
```

