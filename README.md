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
