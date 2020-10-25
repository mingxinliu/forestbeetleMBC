#!/usr/bin/env nextflow

/**
	nfp4MBC - A nextflow pipeline for metabarcoding amplicon sequences
	Copyright (C) 2019	Mingxin Liu - mingxin.liu@utas.edu.au
*/

/**
	Version information
	Version 1.01.10.2019
	Version 2.20.03.2020
*/
	

/**
	Print help when asked for
*/
if (params.help){
	System.out.println("")
	System.out.println("nfp4MBC - A nextflow pipeline for metabarcoding amplicon sequences - Version: 1.20.03.2020 (20.03.2020)")
	System.out.println("")
	System.out.println("Usage: ")
	System.out.println("	nextflow run nfp4MBC.nf ")
	exit 1
}

// Step 0. - Check input parameters

// Checking user-defined parameters
if (params.marker != "CO1" && params.marker != "16S"){
        exit 1, "Marker is not supported! Choose from  <CO1|16S>"
}

if (params.cluster != "usearch" && params.cluster != "swarm"){
	exit 1, "Clustering option not available! Choose from  <usearch|swarm>"
}

// Create working dir
workingpath = params.outdir + "/output_" + params.marker + "_" + params.length
workingdir = file(workingpath)
if (!workingdir.exists()){
	if (!workingdir.mkdirs()){
		exit 1, "Cannot create working directory: $workingpath"
	}
}

// Read paired fastq files
if (params.marker == "CO1"){
	read_pairs = Channel.fromFilePairs("$PWD/fq_${params.marker}/*_L001_R{1,2}_001.fastq", flat: true).ifEmpty("No fastq files are found")
}
if (params.marker == "16S"){
	read_pairs = Channel.fromFilePairs("$PWD/fq_${params.marker}/*_L001_R{1,2}_001.fastq", flat: true).ifEmpty("No fastq files are found")
}
else{
	Channel.empty().ifEmpty("No fastq files are provided")
}

// Step 1.0 - Merge forward and reverse reads with usearch
process merge {
	publishDir "${workingdir}/output_${params.marker}_1.0_merged", mode: "copy", overwrite: true	

	input:
	set sample_id, file (rawForward), file (rawReverse) from read_pairs

	output:
	set sample_id, file ("${sample_id}_merged.fastq") into (merged_c1, merged_c2)

	script:
	"""
	usearch -threads ${params.threads} -fastq_mergepairs $rawForward -fastqout ${sample_id}_merged.fastq
	"""
}

// Step 1.1 - QC with fastqc for merged reads
process fastqc_merged {
	publishDir "${workingdir}/output_${params.marker}_1.1_fastqc_merged", mode: "copy", overwrite: true
	
	input:
	set sample_id, file (in_merged_c1) from merged_c1

	output:
	file ("${sample_id}_fastqc_merged") into fastqc_merged
	
	"""
	mkdir ${sample_id}_fastqc_merged
	fastqc -t ${params.threads} -f fastq -q --casava $in_merged_c1 -o ${sample_id}_fastqc_merged
	"""
}

// Step 1.2- QC with multiqc
process multiqc_merged {
	publishDir "${workingdir}/output_${params.marker}_1.2_multiqc_merged", mode: "copy", overwrite: true
	
	input:
	file ("*") from fastqc_merged.collect().ifEmpty([])

	output:
	file ("multiqc_report.html")
	file "multiqc_data"	

	"""
	multiqc .
	"""
}

// Step 2.0 - remove forward and reverse adapters with CUTADAPT
process cutadapters {
	publishDir "${workingdir}/output_${params.marker}_2.0_cutadatpers", mode: "copy", overwrite: true	

	input:
	set sample_id, file (in_merged_c2) from merged_c2

	output:
	set sample_id, file ("${sample_id}_adapters_trimmed.fastq") into (adapters_trimmed_c1, adapters_trimmed_c2)
	set sample_id, file ("${sample_id}_adapters_untrimmed.fastq")
	
	"""
	cutadapt -a file:$PWD/adapters_${params.marker}/${sample_id}_adapters.fasta -o ${sample_id}_adapters_trimmed.fastq --untrimmed-o ${sample_id}_adapters_untrimmed.fastq $in_merged_c2
	"""
}

// Step 2.1 - report quality control with fastqc
process fastqc_cutadapters {
        publishDir "${workingdir}/output_${params.marker}_2.1_fastqc_cutadapters", mode: "copy", overwrite: true

        input:
        set sample_id, file (in_adapters_trimmed_c1) from adapters_trimmed_c1

        output:
        file ("${sample_id}_fastqc_cutadapters") into fastqc_cutadapters

        """
        mkdir ${sample_id}_fastqc_cutadapters
        fastqc -t ${params.threads} -f fastq -q --casava $in_adapters_trimmed_c1 -o ${sample_id}_fastqc_cutadapters
        """
}

// Step 2.2 - report quality control with multiqc
process multiqc_cutadapters {
        publishDir "${workingdir}/output_${params.marker}_2.2_multiqc_cutadapters", mode: "copy", overwrite: true

        input:
        file ("*") from fastqc_cutadapters.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
	   file "multiqc_data"

        """
        multiqc .
        """
}

// Step 3.0 - remove forward and reverse primers with cutadapt
process cutprimers {
	publishDir "${workingdir}/output_${params.marker}_3.0_cutprimers", mode: "copy", overwirte: true

	input:
	set sample_id, file (in_adapters_trimmed_c2) from adapters_trimmed_c2

	output:
	set sample_id, file ("${sample_id}_primers_trimmed.fastq") into (primers_trimmed_c1, primers_trimmed_c2)
	set sample_id, file ("${sample_id}_primers_untrimmed.fastq")

	script:
	if (params.marker == "CO1")
		"""
		cutadapt --no-indels -M ${params.M} -m ${params.m} -a "AGATATTGGAACWTTATATTTTATTTTTGG;max_error_rate=0.04;min_overlap=29"..."GGAGGATTTGGWAATTGATTAGTW\$;max_error_rate=0.1;min_overlap=22" -o ${sample_id}_primers_trimmed.fastq --untrimmed-o ${sample_id}_primers_untrimmed.fastq $in_adapters_trimmed_c2
		"""

	else if ( params.marker == "16S")
		"""
		cutadapt --no-indels -M ${params.M} -m ${params.m} -a "AGACGAGAAGACCCTATAGA;max_error_rate=0.05;min_overlap=19"..."TACCTTAGGGATAACAGCGTA\$;max_error_rate=0.05;min_overlap=20;" -o ${sample_id}_primers_trimmed.fastq --untrimmed-o ${sample_id}_primers_untrimmed.fastq $in_adapters_trimmed_c2
		"""

	else
		error "Invalid marker!"
}

// Step 3.1 - report quality control with fastqc
process fastqc_cutprimers {
        publishDir "${workingdir}/output_${params.marker}_3.1_fastqc_cutprimers", mode: "copy", overwrite: true

        input:
        set sample_id, file (in_primers_trimmed_c1) from primers_trimmed_c1

        output:
        file ("${sample_id}_fastqc_cutprimers") into fastqc_cutprimers

        """
        mkdir ${sample_id}_fastqc_cutprimers
        fastqc -t ${params.threads} -f fastq -q --casava $in_primers_trimmed_c1 -o ${sample_id}_fastqc_cutprimers
        """
}

// Step 3.2 - report quality control with multiqc
process multiqc_cutprimers {
        publishDir "${workingdir}/output_${params.marker}_3.2_multiqc_cutprimers", mode: "copy", overwrite: true

        input:
        file ("*") from fastqc_cutprimers.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
	file "multiqc_data"

        """
        multiqc .
        """
}

// Step 4.0 - filtering clean reads with usearch
process filter {
	publishDir "${workingdir}/output_${params.marker}_4.0_filtered", mode: "copy", overwrite: true

	input:
	set sample_id, file (in_primers_trimmed_c2) from primers_trimmed_c2

	output:
	file ("${sample_id}_filtered.fasta") into (filtered_fasta_c1, filtered_fasta_c2, filtered_fasta_c3)
	set sample_id, file ("${sample_id}_filtered.fastq") into filtered_fastq_c1
	
	"""
	usearch -threads ${params.threads} -fastq_filter $in_primers_trimmed_c2 -fastq_maxee 1.0 -fastaout ${sample_id}_filtered.fasta -fastqout ${sample_id}_filtered.fastq -relabel @
	"""
}

// Step 4.1 - report quality control with fastqc
process fastqc_filtered {
        publishDir "${workingdir}/output_${params.marker}_4.1_fastqc_filtered", mode: "copy", overwrite: true

        input:
        set sample_id, file (in_filtered_fastq_c1) from filtered_fastq_c1

        output:
        file ("${sample_id}_fastqc_filtered") into fastqc_filtered

        """
        mkdir ${sample_id}_fastqc_filtered
        fastqc -t ${params.threads} -f fastq -q --casava $in_filtered_fastq_c1 -o ${sample_id}_fastqc_filtered
        """
}

// Step 4.2 - report quality control with multiqc
process multiqc_filtered {
        publishDir "${workingdir}/output_${params.marker}_4.2_multiqc_filtered", mode: "copy", overwrite: true

        input:
        file ("*") from fastqc_filtered.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
	file "multiqc_data"

        """
        multiqc .
        """
}

// Step 5.0 - concatenate all filtered reads into one file and find unique reads
process uniques {
	publishDir "${workingdir}/output_${params.marker}_5.0_uniques", mode: "copy", overwrite: true

	input:
	file ("*.filtered.fasta") from filtered_fasta_c1.toList()

	output:
	file ("uniques_${params.marker}.fasta") into (all_uniques)
	
	"""
	cat *.filtered.fasta > all_filtered.fasta
	usearch -fastx_uniques all_filtered.fasta -fastaout uniques_${params.marker}.fasta -relabel Uniq -sizeout -minuniquesize 2 -threads ${params.threads}
	"""
}

// Step 6.0 - cluster all unique reads with user-defined algorithm
process cluster {
	publishDir "${workingdir}/output_${params.marker}_6.0_clustered", mode: "copy", overwrite: true
	
	input:
	file ("uniques_${params.marker}.fasta") from all_uniques

	output:
	file ("zotus_${params.marker}.fasta") into (zotus_rep1, zotus_rep2, zotus_rep3, zotus_rep4, zotus_rep5)

	script:
	if(params.cluster == 'usearch')
	"""
	usearch -unoise3 uniques_${params.marker}.fasta -zotus zotus_${params.marker}.fasta -tabbedout zotus_${params.marker}.txt
	"""

	else if(params.cluster == 'vsearch')
	"""
	vsearch --cluster_unoise uniques_${params.marker}.fasta --centroids zotus_${params.marker}.fasta  --threads ${params.threads}
	"""

	else
	error "Invalid clustring mode: ${mode}"
}


// Step 7.0 - CO1 - map sample filtered reads into clustered zotus with usearch
process zotutab_CO1 {
	publishDir "${workingdir}/output_${params.marker}_7.0_otutab", mode: "copy", overwrite: true

	input:
	file ("*.filtered.fasta") from filtered_fasta_c2.toList()
	file ("zotus_${params.marker}.fasta") from zotus_rep1

	output:
	file ("zotutab_${params.marker}.txt") into zotutab_CO1

	when:
	!params.skip	

	"""
	cat *.filtered.fasta > all_filtered.fasta
	usearch -otutab all_filtered.fasta -sample_delim . -zotus zotus_${params.marker}.fasta -id 0.97 -otutabout zotutab_${params.marker}.txt -threads ${params.threads}
	"""
}

// Step 7.0 - 16S - map sample filtered reads into clustered zotus with usearch
process zotutab_16S {
	publishDir "${workingdir}/output_${params.marker}_7.0_otutab", mode: "copy", overwrite: true

	input:
	file ("*.filtered.fasta") from filtered_fasta_c3.toList()
	file ("zotus_${params.marker}.fasta") from zotus_rep1

	output:
	file ("zotutab_${params.marker}.txt") into zotutab_16S

	when:
	params.skip	

	"""
	cat *.filtered.fasta > all_filtered.fasta
	usearch -otutab all_filtered.fasta -sample_delim . -zotus zotus_${params.marker}.fasta -id 0.97 -otutabout zotutab_${params.marker}.txt -threads ${params.threads}
	"""
}

// Step 8.0 - CO1 - prepare blast matchlist for lulu
process matchlist_CO1 {
	publishDir "${workingdir}/output_${params.marker}_8.0_matchlist", mode: "copy"

	input:
	file ("zotus_${params.marker}.fasta") from zotus_rep2

	output:
	file ("matchlist_${params.marker}.txt") into matchlist_CO1

	when:
	!params.skip

	"""
	makeblastdb -in zotus_${params.marker}.fasta -parse_seqids -dbtype nucl	
	blastn -db zotus_${params.marker}.fasta -outfmt '6 qseqid sseqid pident' -out matchlist_${params.marker}.txt -qcov_hsp_perc 80 -perc_identity 84 -query zotus_${params.marker}.fasta
	"""
}

// Step 8.0 - 16S - prepare blast matchlist for lulu
process matchlist_16S {
	publishDir "${workingdir}/output_${params.marker}_8.0_matchlist", mode: "copy"

	input:
	file ("zotus_${params.marker}.fasta") from zotus_rep2

	output:
	file ("matchlist_${params.marker}.txt") into matchlist_16S

	when:
	params.skip

	"""
	makeblastdb -in zotus_${params.marker}.fasta -parse_seqids -dbtype nucl	
	blastn -db zotus_${params.marker}.fasta -outfmt '6 qseqid sseqid pident' -out matchlist_${params.marker}.txt -qcov_hsp_perc 80 -perc_identity 84 -query zotus_${params.marker}.fasta
	"""
}

// Step 9.0 - CO1 - process zotus with LULU in R
process lulu_CO1 {
	publishDir "${workingdir}/output_${params.marker}_9.0_lulu", mode: "copy", overwrite: true

	input:
	file ("zotutab_${params.marker}.txt") from zotutab_CO1
	file ("matchlist_${params.marker}.txt") from matchlist_CO1

	output:
	file ("zotus_lulu_${params.marker}.txt") into zotus_lulu_CO1
	file ("zotutab_lulu_${params.marker}.csv") into zotutab_lulu_CO1

	"""
	Rscript --vanilla $PWD/lulu_CO1.R zotutab_${params.marker}.txt matchlist_${params.marker}.txt zotus_lulu_${params.marker}.txt zotutab_lulu_${params.marker}.csv
	"""
}

// Step 9.0 - 16S - process zotus with LULU in R
process lulu_16S {
        publishDir "${workingdir}/output_${params.marker}_9.0_lulu", mode: "copy", overwrite: true

        input:
        file ("zotutab_${params.marker}.txt") from zotutab_16S
        file ("matchlist_${params.marker}.txt") from matchlist_16S

        output:
        file ("zotus_lulu_${params.marker}.txt") into zotus_lulu_16S
        file "zotutab_lulu_${params.marker}.csv"

        """
        Rscript --vanilla $PWD/lulu_CO1.R zotutab_${params.marker}.txt matchlist_${params.marker}.txt zotus_lulu_${params.marker}.txt zotutab_lulu_${params.marker}.csv
        """
}

// Step 10.0 - 16S - subseq of lulu curated zotus and blastn against reference sequences
process blastn_16S {
	publishDir "${workingdir}/output_${params.marker}_10.0_blastn", mode: "copy", overwrite: true

	input:
	file ("zotus_lulu_id.txt") from zotus_lulu_16S
	file ("zotus_${params.marker}.fasta") from zotus_rep3
	
	output:
	file ("zotus_lulu_blastn_${params.marker}.txt")	
	
	when:
	params.skip

	"""
	grep "Zotu" zotus_lulu_id.txt | sed 's/^.*Z/Z/' | sed 's/\"//' > lulu_ids.txt
	seqtk subseq zotus_${params.marker}.fasta lulu_ids.txt > zotus_lulu.fasta
	blastn -db $PWD/database/Coleoptera_COI_db -query zotus_lulu.fasta -evalue 1e-4 -outfmt 5 -max_target_seqs 50 -out zotus_lulu_blastn_${params.marker}.txt -num_threads ${params.threads}
	"""
}

// Step 10.1 - CO1 - remove zotus with inframe stop codons
process translate_CO1 {
	publishDir "${workingdir}/output_${params.marker}_10.1_translate", mode: "copy", overwrite: true
	frame = ['1', '2', '3']	

	input:
	file ("zotus_${params.marker}.fasta") from zotus_rep3
	each codon from frame	
	
	output:
	file ("${codon}.txt") into header
	
	when:
	!params.skip	

	"""
	seqkit translate -T 5 -f ${codon} zotus_${params.marker}.fasta > zotus_${params.marker}_AA${codon}.fasta
	grep -B1 \\* zotus_${params.marker}_AA${codon}.fasta | grep -vFf - zotus_${params.marker}_AA${codon}.fasta | sed 's/_frame=${codon}//' | awk 'sub(/^>/, "")' > ${codon}.txt
	"""
}

// Step 10.2 - CO1 - subset sequences with no stop codons
process nostopcodon_CO1 {
	publishDir "${workingdir}/output_${params.marker}_10.2_nostopcodon", mode: "copy", overwrite: true

	input:
	file ("zotus_${params.marker}.fasta") from zotus_rep4
	file ("*.txt") from header.toList()
		
	output:
	file ("zotus_${params.marker}_codon.fasta") into zotus_codon

	when:
	!params.skip

	"""
	cat *.txt > names.txt
	seqtk subseq zotus_${params.marker}.fasta names.txt > subseq.fasta		
	seqkit rmdup -n subseq.fasta -o zotus_${params.marker}_codon.fasta 
	"""
}

// Step 11.0 - CO1 - find common zotus between nostopcodon zotus and lulu curated zotus
process commonzotus_CO1 {
	publishDir "${workingdir}/output_${params.marker}_11.0_commonzotus", mode: "copy", overwrite: true

	input:
	file ("zotus_${params.marker}.fasta") from zotus_rep5
	file ("zotus_lulu_${params.marker}.txt") from zotus_lulu_CO1
	file ("zotus_${params.marker}_codon.fasta") from zotus_codon

	output:
	file ("zotus_lulu_${params.marker}.fasta") into zotus_lulu
	file ("zotus_common_${params.marker}.fasta") into zotus_common
	file ("zotus_common_blastn_${params.marker}.txt")

	when:
	!params.skip

	"""
	grep "Zotu" zotus_lulu_${params.marker}.txt | sed 's/^.*Z/Z/' | sed 's/\"//' > lulu_ids.txt
	seqtk subseq zotus_${params.marker}.fasta lulu_ids.txt > zotus_lulu_${params.marker}.fasta
	seqkit common zotus_lulu_${params.marker}.fasta zotus_${params.marker}_codon.fasta -o zotus_common_${params.marker}.fasta
	blastn -db $PWD/database/Coleoptera_COI_db -query zotus_common_${params.marker}.fasta -evalue 1e-4 -outfmt 5 -max_target_seqs 50 -out zotus_common_blastn_${params.marker}.txt -num_threads ${params.threads}
	"""
}
