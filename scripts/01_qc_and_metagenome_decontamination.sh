#!/bin/bash

#QC#
#Trimmomatic v0.33#
#Perform on both metagenomic and metatranscriptomic reads#
java -jar trimmomatic-0.33.jar PE -trimlog trimmomatic_log_file.txt \
input_forward_read.fastq.gz input_reverse_read.fastq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE-modified.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Decontamination#
#Co-assembly of negative control reads#
#Megahit v1.2.9#
#Use zcat to unzip reads to prepare for megahit
megahit -1 output_forward_neg_sampl_1_paired.fq,output_forward_neg_sampl_2_paired.fq \
-2 output_reverse_neg_sampl_1_paired.fq,output_reverse_neg_sampl_2_paired.fq \
-o /output/folder

#Removal of reads deemed to be contaminants#
#First need to create a database of negative control reads which have already been coassembled#
perl deconseq.pl -f coassembled_negative_controls.fq -dbs decontam_database

#Run DeconSeq#
perl deconseq.pl -f output_forward_paired.fq -dbs decontam_database -c 95 -out_dir /output_directory
perl deconseq.pl -f output_reverse_paired.fq -dbs decontam_database -c 95 -out_dir /output_directory
#-c is the identity cutoff for a match
