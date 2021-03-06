# Bioinfo-Scripts
Record of scripts developed and used by the Kovach lab for processing sequencing data.

Thus far, the repository includes:

1.) BioinformaticBookkeepingSlurmTemplate.sh - A sample slurm scheduler template for generating a bioinformatic "paper trail" and renaming slurm output records more informatively.

2.) DataPreProcessingScript.sh - A script for creating fastqc reports on newly sequenced data, optional code for double checking the adaptors used if needed, and for quality trimming adaptors and collapsing overlapping read pairs with AdapterRemoval2 on UNH's Premise cluster.

3.) BWA-MEM_Alignment.sh - Runs pipeline for mapping quality-trimmed fastq reads to a reference genome with the BWA-MEM algorithm and performs additional quality control steps based on the GATK Best Practices Data Pre-Processing for Variant Discovery Pipeline (including samtools fixmate, index & sort; GATK Indel Realigner, Picard Mark Duplicates, and GATK Base Quality Score Recalibration) to generate analysis-ready bam files. 

4.) AdapterRemovalSettingsCharts.r - an R-based subscript called from DataPreProcessingScript.sh (from version 1.2.2+) that generates nice pdf histogram line charts of the read length distribution data generated by AdapterRemoval.

5.) LEF-FinalDataPreProcessingSlurmScript.sh - A copy of my slurm job scheduler for running the DataPreProcessingScript.sh. (Has additional functionality than the basic BioinformaticBookkeepingSlurmTemplate.sh).
