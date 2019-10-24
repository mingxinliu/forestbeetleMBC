#!/usr/bin/env nextflow

/**
	A nextflow pipeline for processing DNA metabarcoding amplicon sequences
	Copyright (C) 2019	Mingxin Liu - mingxin.liu@utas.edu.au
*/

/**
	Print version when asked for
*/

version='1.0.10.01'
timestamp='20191001'

if (params.version) {
	System.out.println("")
	System.out.println("Mingxin's nextflow pipeline for processing DNA metabarcoding amplicon sequences - version: $version ($timestamp)")
	exit 1
}

/**
	Print help info when asked for
*/

Channel
	.fromFilePairs( "fq_CO1/*_L001_R{1,2}_001.fastq", flat:true )
	.ifEmpty { error "Cannot find any matching: ${params.reads}" }
	.set { read_pairs }

// Step 1.0 - Merge forward and reverse reads with USEARCH
process merge {
	publishDir './output_CO1_1.0_merged', mode: 'copy', overwrite: true

	input:
	set sample_id, file (rawForward), file(rawReverse) from read_pairs

	output:
	set sample_id, file ("${sample_id}.merged.fastq") into (merged_reads_1, merged_reads_2)

	"""
	usearch -threads 12 -fastq_mergepairs $rawForward -fastqout ${sample_id}.merged.fastq
	"""
}

// Step 1.1 - report quality control with fastqc for merged reads
process fastqc_merged {
	publishDir './output_CO1_1.1_fastqc_merged', mode: 'copy', overwrite: true
	
	input:
	set sample_id, file (merged_reads_1) from merged_reads_1

	output:
	file ("${sample_id}_fastqc_merged") into fastqc_merged
	
	"""
	mkdir ${sample_id}_fastqc_merged
	fastqc -threads 12 -o ${sample_id}_fastqc_merged -f fastq -q $merged_reads_1
	"""
}

// Step 1.2- report quality control with multiqc
process multiqc_merged {
	publishDir './output_CO1_1.2_multiqc_merged', mode: 'copy', overwrite: true
	
	input:
	file ("*") from fastqc_merged.collect().ifEmpty([])

	output:
	file ("multiqc_report.html")
	file "multiqc_data" into multiqc_merged

	"""
	multiqc .
	"""
}

// Step 2.0 - remove forward and reverse adapters with CUTADAPT
process cutadapters {
	publishDir './output_CO1_2.0_adapters_trimmed', mode: 'copy', overwrite: true

	input:
	set sample_id, file (merged_reads_2) from merged_reads_2

	output:
	set sample_id, file ("${sample_id}.adapters_trimmed.fastq") into (adapters_trimmed_1, adapters_trimmed_2)
	
	"""
	cutadapt -a file:$PWD/adapters_CO1/${sample_id}_adapters.fasta -o ${sample_id}.adapters_trimmed.fastq $merged_reads_2
	"""
}

// Step 2.1 - report quality control with fastqc
process fastqc_cutadapters {
        publishDir './output_CO1_2.1_fastqc_cutadapters', mode: 'copy', overwrite: true

        input:
        set sample_id, file (adapters_trimmed_1) from adapters_trimmed_1

        output:
        file ("${sample_id}_fastqc_cutadapters") into fastqc_cutadapters

        """
        mkdir ${sample_id}_fastqc_cutadapters
        fastqc -threads 12 -o ${sample_id}_fastqc_cutadapters -f fastq -q $adapters_trimmed_1
        """
}

// Step 2.2 - report quality control with multiqc
process multiqc_cutadapters {
        publishDir './output_CO1_2.2_multiqc_cutadapters', mode: 'copy', overwrite: true

        input:
        file ("*") from fastqc_cutadapters.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
        file "multiqc_data" into multiqc_cutadapters

        """
        multiqc .
        """
}

// Step 3.0 - remove forward and reverse primers with CUTADAPT
process cutprimers {
	publishDir './output_CO1_3.0_primers_trimmed', mode: 'copy', overwrite: true

	input:
	set sample_id, file (adapters_trimmed_2) from adapters_trimmed_2

	output:
	set sample_id, file ("${sample_id}.primers_trimmed.fastq") into (primers_trimmed_1, primers_trimmed_2)

	"""
	cutadapt --no-indels -M 157 -m 157 -a "AGATATTGGAACWTTATATTTTATTTTTGG;max_error_rate=0.04;min_overlap=29"..."GGAGGATTTGGWAATTGATTAGTW\$;max_error_rate=0.1;min_overlap=22" -o ${sample_id}.primers_trimmed.fastq $adapters_trimmed_2
	"""
}

// Step 3.1 - report quality control with fastqc
process fastqc_cutprimers {
        publishDir './output_CO1_3.1_fastqc_cutprimers', mode: 'copy', overwrite: true

        input:
        set sample_id, file (primers_trimmed_1) from primers_trimmed_1

        output:
        file ("${sample_id}_fastqc_cutprimers") into fastqc_cutprimers

        """
        mkdir ${sample_id}_fastqc_cutprimers
        fastqc -threads 12 -o ${sample_id}_fastqc_cutprimers -f fastq -q $primers_trimmed_1
        """
}

// Step 3.2 - report quality control with multiqc
process multiqc_cutprimers {
        publishDir './output_CO1_3.2_multiqc_cutprimers', mode: 'copy', overwrite: true

        input:
        file ("*") from fastqc_cutprimers.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
        file "multiqc_data" into multiqc_cutprimers

        """
        multiqc .
        """
}

// Step 4.0 - filtering clean reads with USEARCH
process filter {
	publishDir './output_CO1_4.0_filtered', mode: 'copy', overwrite: true

	input:
	set sample_id, file (primers_trimmed_2) from primers_trimmed_2

	output:
	file ("${sample_id}.filtered.fasta") into (filtered_reads_1, filtered_reads_2)
	set sample_id, file ("${sample_id}.filtered.fastq") into filtered_reads_3
	
	"""
	usearch -threads 12 -fastq_filter $primers_trimmed_2 -fastq_maxee 1.0 -fastaout ${sample_id}.filtered.fasta -fastqout ${sample_id}.filtered.fastq -relabel @
	"""
}

// Step 4.1 - report quality control with fastqc
process fastqc_filtered {
        publishDir './output_CO1_4.1_fastqc_filtered', mode: 'copy', overwrite: true

        input:
        set sample_id, file (filtered_reads_3) from filtered_reads_3

        output:
        file ("${sample_id}_fastqc_filtered") into fastqc_filtered

        """
        mkdir ${sample_id}_fastqc_filtered
        fastqc -threads 12 -o ${sample_id}_fastqc_filtered -f fastq -q $filtered_reads_3
        """
}

// Step 4.2 - report quality control with multiqc
process multiqc_filtered {
        publishDir './output_CO1_4.2_multiqc_filtered', mode: 'copy', overwrite: true

        input:
        file ("*") from fastqc_filtered.collect().ifEmpty([])

        output:
        file ("multiqc_report.html")
        file "multiqc_data" into multiqc_filtered

        """
        multiqc .
        """
}

// Step 5.0 - concatenate all filtered reads into one file and find unique reads
process uniques {
	publishDir './output_CO1_5.0_uniques', mode: 'copy', overwrite: true

	input:
	file ("*.filtered.fasta") from filtered_reads_1.toList()

	output:
	file ("all_uniques.fasta") into (unique_reads_A)
	
	"""
	cat *.filtered.fasta > all_filtered_1.fasta
	usearch -fastx_uniques all_filtered_1.fasta -fastaout all_uniques.fasta -relabel Uniq -sizeout -minuniquesize 2 -threads 12
	"""
}

// Step 6.0 - cluster all unique reads with USEARCH
process cluster_USEARCH {
	publishDir './output_CO1_6.0_otus_USEARCH', mode: 'copy', overwrite: true
	
	input:
	file ("all_uniques.fasta") from unique_reads_A

	output:
	file ("otus_USEARCH.fasta") into (otus_USEARCH_1, otus_USEARCH_2)
	file	"unoise3.txt"	
	
	"""
	usearch -unoise3 all_uniques.fasta -zotus otus_USEARCH.fasta -tabbedout unoise3.txt
	"""
}

// Step 7.0 - map sample filtered reads into clustered zotus with USEARCH
process otutab {
	publishDir './output_CO1_7.0_otutab', mode: 'copy', overwrite: true

	input:
	file ('*.filtered.fasta') from filtered_reads_2.toList()
	file ("otus_USEARCH.fasta") from otus_USEARCH_1

	output:
	file ("otutab.txt") into otutab

	"""
	cat *.filtered.fasta > all_filtered_2.fasta
	usearch -otutab all_filtered_2.fasta -sample_delim . -zotus otus_USEARCH.fasta -id 0.97 -otutabout otutab.txt -threads 12
	"""
}


// Step 8.0 - prepare blast matchlist for LULU
process matchlist {
	publishDir './output_CO1_8.0_LULU', mode: 'copy', overwrite: true

	input:
	file ("otus_USEARCH.fasta") from otus_USEARCH_1

	output:
	file ("matchlist.txt") into matchlist

	"""
	makeblastdb -in otus_USEARCH.fasta -parse_seqids -dbtype nucl	
	blastn -db otus_USEARCH.fasta -outfmt '6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc 80 -perc_identity 84 -query otus_USEARCH.fasta
	"""
}

// Step 8.1 - process zotus with LULU in R
process lulu {
	publishDir './output_CO1_8.0_LULU', mode: 'copy'

	input:
	file ("otutab.txt") from otutab	
	file ("matchlist.txt") from matchlist
	
	output:
	file ("lulu_otus_id.txt") into lulu_otus_id
	file ("lulu_result.csv")

	"""
	Rscript --vanilla $PWD/lulu_CO1.R otutab.txt matchlist.txt lulu_otus_id.txt lulu_result.csv
	"""
}

// Step 9.0 - BLAST against reference sequence database and output file for MEGAN
process blastn {
	publishDir './output_CO1_9.0_MEGAN', mode: 'copy', overwrite: true

	input:
	file ("lulu_otus_id.txt") from lulu_otus_id
	file ("otus_USEARCH.fasta") from otus_USEARCH_2

	output:
	file "lulu_ids.txt"
	file "lulu_otus.fasta"
	file "lulu_otus_MEGAN.txt"
	
	"""
	grep "Zotu" lulu_otus_id.txt | sed 's/^.*Z/Z/' | sed 's/\"//' > lulu_ids.txt
	seqtk subseq otus_USEARCH.fasta lulu_ids.txt > lulu_otus.fasta
	blastn -db $PWD/database_CO1/Coleoptera_CO1.fasta -query lulu_otus.fasta -evalue 0.001 -outfmt 5 -max_target_seqs 50 -out lulu_otus_MEGAN.txt -num_threads 12
	"""
}
