#!/usr/bin/env bash

# Author: Ana Herráez Pérez
# Date : October 2026
# Language : Bash
# Description : Semi-automated assembly-processing pipeline designed for Africonus crotchii. Input : PacBio CCS (HiFi) fastq. Output : assembly.
# Publication : Herráez-Pérez, A., Pardos-Blas, J.R., Zardoya, R. 2026. Haplotype-Resolved Chromosome-Level Genome of the Cabo Verde Venomous Snail Africonus crotchii. G3

###############################
### De novo genome assembly ###
###############################

# Hifiasm v0.15.4-r347 - Haplotype-resolved de novo assembly
	hifiasm -o hifiasm.asm -t32 --h1 hifi_read1.fq.gz --h2 hifi_read2.fq.gz 2> hifiasm.log

# Blast v2.13.0 - Detection of contamination
	blastn \
	-query primary.fa \
	-db nt \
	-outfmt '6 qseqid staxids bitscore std' \
	-max_target_seqs 10 \
	-max_hsps 1 \
	-evalue 1e-25 \
	-num_threads 32 \
	-out blast_nt.out
	# Independent primary.fa for haplotype one and two

# Minimap2 v2.21 - Read mapping
	minimap2 -ax map-hifi -t 32 primary.fa reads.fastq.gz | samtools sort -@ 16 -o reads.sorted.bam
	samtools index reads.sorted.bam
	
# blobtools v1.1.1 - Taxonomic assessment
	blobtools create -i primary.fa -b reads.sorted.bam -t blast_nt.out -o blob_out
	blobtools view -i blob_out.blobDB.json -r phylum superkingdom -o blob_out
	blobtools plot -i blob_out.blobDB.json -o blob_out
	
# Seqkit - Remove contaminant reads
	seqkit grep -v -f contaminants.txt primary.fa > primary.filtered.fa

#######################################################
## Genome scaffolding using proximity ligation data ###
#######################################################

#	BWA -  Omni-C alignment
	# Inputs
	hap1=hap1.purged.fa
	hap2=hap2.purged.fa

	R1=reads_R1.fastq.gz
	R2=reads_R2.fastq.gz

	# Index both assemblies
	bwa index "$hap1"
	bwa index "$hap2"

	# Haplotype one
	bwa mem -5SP -T0 -t 32 "$hap1" "$R1" "$R2" | samtools sort -@ 16 -o hap1.sorted.bam
	samtools index hap1.sorted.bam	

	# Haplotype 2
	bwa mem -5SP -T0 -t 32 "$hap2" "$R1" "$R2" | samtools sort -@ 16 -o hap2.sorted.bam
	samtools index hap2.sorted.bam
	
	# Count primary alignments with MAPQ > 50
	samtools view -c -q 51 -F 0x904 hap1.sorted.bam
	samtools view -c -q 51 -F 0x904 hap2.sorted.bam
	
# HiRise - Software pipeline designed specifi cally for usingproximity ligation data to scaffold genome assemblies
	# Performed independently for each haplotype using Omni-C reads.
	# Only read pairs with mapping quality >50 were considered.
	# No command line

# Converting fastq to fasta
	sed -n '1~4s/^@/>/p;2~4p' input.fastq.gz > output.fasta

# Flye - raw assembly
	flye --nano-raw output.filtlong.fasta --genome-size 23m --threads 20 --asm-coverage 35 --debug --out-dir output_flye 
	#Threshold to be adjusted according to sequencing results: "--asm-coverage" 

# Quast - assembly quality check
	quast input_contig.fasta -r input-reference.fasta --threads 10 --labels "title" -o output_quast --glimmer --nanopore input.filtlong.fastq.gz 

########################
## Genome annotation ###
########################

#1# Repeat annotation - All steps wereperformed indepently for each haplotype

# RepeatModeler v2.0.6 - De novo transposable element (TE) family identification and modeling
	# Build database
	BuildDatabase -name TE_DB hic_output.fasta
	# TE identification
	RepeatModeler -threads 15 -database TE_DB -LTRStruct

# RepeatMasker v4.1.7. - Screens DNA sequences for interspersed repeats and low complexity DNA sequences
	LAST_DIR=$(ls -dt */ | head -n 1) # RepeatModeler output directory
	RepeatMasker hic_output.fasta -pa 25 -gff \
	-xsmall -s -gccalc -a \
	-engine ncbi \
	-lib $LAST_DIR/consensi.fa.classified \
	-dir RMaskerOutput

#2# RNA-processing

# Trimmomatic - Clean RNA-seq reads
	java -jar Trimmomatic-0.39.jar PE \
	reads_R1.fastq.gz reads_R2.fastq.gz \
	trimmed_R1.fastq.gz unpaired_R1.fastq.gz \
	trimmed_R2.fastq.gz unpaired_R2.fastq.gz \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
	#Repeat this command for all samples
	
# STAR v2.7.11b - RNA-seq aligment to a reference
	# All steps were performed independently for each haplotype.
	# Generate genome index
	mkdir -p star_index
	STAR --runMode genomeGenerate \
	--runThreadN 32 \
	--genomeDir star_index \
	--genomeFastaFiles assembly_softmasked.fa

	# Align trimmed paired-end RNA-seq reads
	STAR --genomeDir star_index \
	--runThreadN 32 \
	--readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix sample_
	# Repeat this command for each RNA-seq sample, using the same genome index for all samples of a given haplotype.
	# Assign a unique output prefix to each sample.

	# Merge RNA-seq alignments with samtools
	# Perform independently for each haplotype.
	# Combine the coordinate-sorted BAM files generated from all RNA-seq samples mapped to the same haplotype.
	# Do not combine BAM files aligned to different haplotypes.
	# Include all available RNA-seq samples in the merge command (sample1_, sample2_, ...).
	samtools merge -@ 16 -f rnaseq.merged.bam \
	sample1_Aligned.sortedByCoord.out.bam \
	sample2_Aligned.sortedByCoord.out.bam \
	sample3_Aligned.sortedByCoord.out.bam

	# Sort and index the merged BAM file.
	samtools sort -@ 16 -o rnaseq.sorted.bam rnaseq.merged.bam
	samtools index rnaseq_sorted.bam

#3# Gene annotation
	# Prepare external protein evidence. Combine the protein datasets from the two related species and Metazoa OrthoDB v12.
	# Ensure that protein sequence identifiers are unique.
	# This step only needs to be performed once. The same protein database can be used for both haplotypes.
	cat protein_species1.fa protein_species2.fa metazoa_orthodb12.fa > proteins_merged.fa
	
#BRAKER v3.0.8 - Uses genomic and RNA-Seq data to automatically generate full gene structure annotations in novel genome.
	# Perform gene prediction independently for each haplotype.
	GE=assembly_softmasked.fa
	BA=rnaseq_sorted.bam
	PR=seq proteins_merged.fa
	
	braker.pl --genome=$GE --bam=$BA --prot_seq=$PR --softmasking --etpmode --cores 16 --species=species

	##Extraction of supported genes
	file1=augustus.hints.gff3
	file2=hintsfile.gff
	
	python selectSupportedSubsets.py --fullSupport fullsupport.gff --anySupport anySupport.gff --noSupport nosupport.gff $file1 $file2
	gffread -x augustus_genes_anysupport.fasta -g $GE anySupport.gff

###############################################
### Quality assessment of genome assemblies ###
###############################################

# All steps were performed independently for each haplotype.

# QUAST v5.0.2. - Genome assembly QUality ASsessment 
	python3 quast.py assembly.fasta -t 10 --large -o quast_results
	
# Minimap2 v2-2.30-r1287 - Map long reads to assembly
	#Variables
	REF=assembly.fasta
	fastq1=hifi_read1.fq.gz
	fastq2=hifi_read2.fq.gz
	T=Threads #Number of threads
	OUT=minimap2_output.bam
	
	# Index genome assembly fasta file
	samtools faidx $REF
	
	# Mapping reads to reference. Output: Sorted and indexed BAM file.
	minimap2 -t $T -ax map-pb $REF $fastq1 $fastq2 | samtools view -Sb - | samtools sort -o $OUT -@20 -
	
	# Calculate coverage depth
	samtools depth -a $OUT > coverage_per_base.txt
	
	#Calculate genome coverage percentage
	GENOME_LENGTH=$(awk '{sum+=$2} END {print sum}' ${REFERENCE}.fai)  # Tamaño total del genoma
	echo "Genome lenght is $GENOME_LENGTH" 
	
	COVERED_BASES=$(awk '{ if($3 > 0) covered++ } END { print covered }' coverage_per_base.txt)
	COVERAGE_PERCENTAGE=$(awk -v covered=$COVERED_BASES -v total=$GENOME_LENGTH 'BEGIN { print (covered / total) * 100 }')
	echo "Genome coverage percentage: $COVERAGE_PERCENTAGE %"

# Mosdepth v0.3.13 - Fast BAM/CRAM depth calculation for WGS
	BAM=minimap2_output.bam
	mosdepth sample $BAM
	
# blobtoolkit v4.1.5 - Quality and contamination assessment of genome assemblies
	# BLASTn
	# Variables
	REF="/path/to/assembly.fasta"
	T=Threads # Set number of threads
	# Path to BLASTN Databases.
	nt="path/to/Archaea-Bacteria-Virus_NCBI.fasta"
	F="path/to/Fungi_genome.fna"
	M="path/to/mollusc_genome_LV_KC.fna" #M=mollusca genome database + Kcan y Lven genomes
	P="path/to/platelminte.fna"
	N="path/to/nematode.fna"
	H="path/to/homos.fna"
	PO="path/to/polichaeta.fna"
	# Path to BAM file (from minimap2)
	BAM="path/to/minimap2_output.bam"
	
	#Running blastn indepently for each database
	blastn -query $REF -db $nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out nt.tsv
	blastn -query $REF -db $F -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out fungi.tsv
	blastn -query $REF -db $M -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out mollusca.tsv
	blastn -query $REF -db $P -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out platelminte.tsv
	blastn -query $REF -db $N -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out nematode.tsv
	blastn -query $REF -db $H -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out homos.tsv
	blastn -query $REF -db $PO -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads $T -out polichaeta.tsv
	
	# Blastn output variables
	nt=nt.tsv
	nematode=fungi.tsv
	platelminte=mollusca.tsv
	fungi=platelminte.tsv
	mollusca=nematode.tsv
	human=homos.tsv
	polychaeta=polichaeta.tsv
	
	# BUSCO output variables
	metazoa=Metazoa_full_table.tsv
	eukaryota=Eukaryota_full_table.tsv
	mollusca=Mollusca_full_table.tsv 
	
	# Blobtools output dir
	OUT=blobtools_out
	# Running blobtools
	blobtools create \
		--fasta $REF \
		--cov $BAM \
		--taxdump /path/to/taxdump \  #Set path directory
		--hits $nt \
		--hits $nematode \
		--hits $platelminte \
		--hits $fungi \
		--hits $mollusca \
		--hits $human \
		--hits $polychaeta \
		--taxrule bestdistsum \
		--busco $metazoa \
		--busco $eukaryota \
		--busco $mollusca \
		$OUT
	
	# Results visualization
	# Let’s suppose that you have obtained a results directory. Now, copy it to you local computer (/path/to/dest_dir) and run this command:
	blobtools view --local /path/to/dest_dir
	# !! Note that this command requires a path to the directory where results dir resides.
	
# Meryl v1.4.1. - Genomic k-mer count
	# Variables
	REF=assembly.fasta
	R1=hifi_read1.fastq.gz
	R2=hifi_read2.fastq.gz
	T=Threads # Set the number of threads
	mem=memory # Set memory. i.e. 128g.
	OUT=meryl_out

	# K-mer counting
	meryl count k=21 threads=$T memory=$mem $OUT read_db.meryl $R1 $R2
	
	# Meryl statistics
	meryl statistics read_db.meryl > read_db.hist

# Merqury v1.3. - Evaluate genome assemblies with k-mers
	# Variables
	H1=hap1_assembly.fasta
	H2=hap2_assembly.fasta
	OUT=merqury_out
	
	# Running merqury
	merqury.sh read_db.meryl $H1 $H2 $OUT

	# Generate spectra plot
	Rscript plot_spectra_cn.R -f merqury.spectra-cn.hist -o merqury.spectra-cn
	
# BUSCO v6.0.0. - Genome assembly and annotation completeness
	# BUSCO assessment was run indepently using metazoa_odb12 and mollusca_odb12
	# Variables
	A=hic_output.fasta #Genome assembly
	P=braker.aa # Protein sequences from genome anotation in fasta file
	B=lineage_odb12 # Set database
	C=N #Set number of threads
	O=busco_out
	config="/path/to/config.ini"
	
	# Gemone mode
	python3 busco --ofline -m genome -i $A -o $O -c $C -l $B --config $config

	# Protein mode (annotation)
	python3 busco --offline -m protein -i $P -o $O -c $C -l $B --config $config


############################
### Haplotype comparison ###
############################

# Gene structure comparison and boxplots
# Gene, exon and intron lengths were extracted from BRAKER3 annotations using a custom Python script: boxplots_genes.py
# The longest transcript was selected for each gene.
# Both haplotypes were compared in the same analysis.

	python3 boxplots_genes.py --hap1 hap1.braker.gtf --hap2 hap2.braker.gtf --outdir gene_structure

# Repeat landscapes
# Run independently for each haplotype.
# Use the corresponding RepeatMasker .align file.
	GENOME_SIZE=$(awk '{sum+=$2} END {print sum}' ${REFERENCE}.fai)  # Tamaño total del genoma

	calcDivergenceFromAlign.pl -s assembly.divsum assembly.fa.align
	createRepeatLandscape.pl -div assembly.divsum -g GENOME_SIZE > assembly.landscape.html

# Whole-genome alignment
	# Haplotype one: query
	# Haplotype two: target
	minimap2 -t 20 -c -x asm5 --cs hap2.fa hap1.fa > hap1_to_hap2.paf

	# Visualization with D-GENIES v1.5.0.
	# Upload the PAF file and the corresponding assemblies.
	# Online platform. No command line

# BUSCO ortholog comparison
	# BUSCO status transitions and genomic coverage of missing orthologs were analyzed using an in-house Python script
	# Examine Missing-to-Complete BUSCOs using primary whole-genome alignment blocks. 
	# Calculate coverage of BUSCO regions and the surrounding +/- 50 kb intervals.
	
	python3 busco_haplotype_comparison.py \
		--hap1-busco hap1/full_table.tsv \
		--hap2-busco hap2/full_table.tsv \
		--paf hap1_to_hap2.paf \
		--hap2-fai hap2.fa.fai \
		--chrom-map chromosome_pairs.tsv \
		--flank 50000 \
		--outdir busco_comparison

########################
### Synteny Analyses ###
########################

# MCScanX v1.0.0. - Analyze syntenic relationships between A. crotchii, K. canariensis, and L. ventricosus
	# Merge protein-coding genes from all four genome assemblies
	cat HAP1_proteins.aa  HAP2_proteins.aa  KCAN_proteins.fasta  LVEN_proteins.fa > all_proteins.fasta
	
	# Make database for diamond
	diamond makedb --in all_proteins.fasta -d all_proteins.dmnd
	
	# All-versus-all searching
	diamond blastp -q all_proteins.fasta -d all_proteins.dmnd -p 5 -e 1e-10 --max-hsps 5 -o all_blastp.tsv
	
	# MCScanX - Multiple Collinearity Scan
	MCScanX all_blastp.tsv
	
	# Visualization using online interactive platform SynVisio (no command line)
	# Upload collinearity and gff files from MCScanX
	
#####################
### End of script ###
#####################
