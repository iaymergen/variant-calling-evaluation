# Variant Calling Pipeline Evaluation

This project evaluates the performance of different variant calling pipelines applied to a germline-tumor paired dataset. The analysis compares combinations of mappers, variant callers, and preprocessing conditions based on metrics like precision, recall, and F1-score.

## 🔧 Pipeline Setup

Twelve pipelines were created by combining:
- 2 mappers: BWA, Bowtie  
- 3 variant callers: Mutect, Strelka, SomaticSniper  
- With/without base recalibration

## 📁 Project Structure
- data/ # FASTQ, reference files, VCF outputs
- scripts/ # Bash and Python scripts
- report_and_results/ # Visualizations, tables, final report


## 📊 Analysis Highlights

- Per-base QC and GC content analysis (FastQC)  
- BAM alignment stats (flagstat, coverage)  
- VCF-level comparison of precision, recall, F1-score  
- Visualizations: PCA, heatmaps, Venn diagrams, runtime comparisons

## 📝 Report
- All findings and visualizations are summarized in report_and_results/Project_Report.pdf

Note: This repository is not designed for direct re-execution.
Some file paths in scripts are project-specific, and large input files (e.g., reference genomes, FASTQ files) are excluded due to size constraints.  
The purpose of this repository is to document the workflow, code structure, and analysis results of the project.

