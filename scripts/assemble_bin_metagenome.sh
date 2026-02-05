#!/bin/bash

###################################################################################################################
#Assemble decontaminated metagenome reads#
#MegaHit v1.2.9#
megahit -1 decontaminated_forward_paired.fq -2 decontaminated_reverse_paired.fq -o /output/folder

###################################################################################################################
#Map metagenome reads to metagenome assembly in prep for binning#
#Bowtie2 v2.5.1#
#build index of metagenome first
bowtie2-build input_assembly.fa output_assembly_index_metagenome
#last bit is setting the first part of the file name the program will use for all the files

#Next actually align the metagenome reads to the indexed metagenome
bowtie2 -x input_assembly_index_metagenome -1 input_decontaminated_forward_reads.fq \
-2 input_decontaminated_reverse_reads.fq -S output_assembly_mapped.sam \
--very-sensitive-local --threads 60
#The sam file is the file that the index function made

#Convert .sam file to .bam file and sort it#
#samtools v1.16.1#
samtools view -bS output_assembly_mapped.sam > output_assembly_mapped.bam
samtools sort output_assembly_mapped.bam -o output_assembly_mapped_sorted.bam

###################################################################################################################
#Bin the metagenome#
#Binning with MetaBat2 v2.15#
runMetaBat.sh -m 2500 assembly_contigs.fa output_assembly_mapped_sorted.bam
#-m sets the minimum contig length.

#Binning with MaxBin v.2.2.7#
perl /home/whytelab_2/miniconda3/opt/MaxBin-2.2.7/run_MaxBin.pl -contig assembly_contigs.fa \
-reads decontaminated_forward_paired.fq \
-reads2 decontaminated_reverse_paired.fq \
-out output_prefix -thread 56

#Binning with SemiBin2 v.1.5.1#
SemiBin2 single_easy_bin --input-fasta assembly_contigs.fa \
--input-bam output_assembly_mapped_sorted.bam --environment global \
--output output_prefix

###################################################################################################################
#Check bin completeness and contamination#
#Perform for bins created using all three methods above#
#CheckM2 v.1.0.0#
checkm2 predict -x fa --input /path/to/input/bins --output-directory /path/to/output/directory

###################################################################################################################
#Taxonomically identify bins#
#GTDB-tk v.2.2.4#
gtdbtk classify_wf -x fa --genome_dir /path/to/input/bins --out_dir /path/to/output/directory

###################################################################################################################
#Dereplication with dRep v.3.4.3#
#First off. dRep uses CheckM to estimate completeness and contamination however,
#I've already used CheckM2 to estimate these values and results can vary a lot between CheckM and CheckM2.
#Luckily, dRep allows you to provide your own completeness and contamination file so using the file guidelines here 
#(https://drep.readthedocs.io/en/master/advanced_use.html#using-external-genome-quality-information),
#I created my own file using my CheckM2 output to input into dRep using the --genomeInfo tag

#To use dRep
dRep dereplicate /path/to/output/directory -g /path/to/bins/*.fa \
--genomeInfo checkM2_results.csv --processors 56

##################################################################################################################
#Estimate in situ bin replication rate#
#First must use bowtie2 to make an index of your bin. Then map metagenome forward and reverse reads to the bin.
#iRep can only handle bins more than 75% complete, less than 2% contaminated and less than 175 scaffold per Mbp.

#Bowtie2 v2.5.1#
#build index of metagenome first
bowtie2-build input_bin.fa output_index_bin
#last bit is setting the first part of the file name the program will use for all the files

#Next actually align the metagenome reads to the indexed bin
bowtie2 -x output_index_bin -1 input_decontaminated_forward_reads.fq \
-2 input_decontaminated_reverse_reads.fq -S output_bin_mapped.sam \
--very-sensitive-local --threads 60

#iRep v.1.10#
iRep -f input_bin.fa -s output_bin_mapped.sam -o iRep_bin.iRep
