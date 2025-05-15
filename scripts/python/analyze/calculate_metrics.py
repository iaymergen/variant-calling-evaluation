import os
import vcfpy  # Use vcfpy instead of vcf
import pandas as pd

# Define paths to ground truth files
vcf_directory = "/home/iaymergen/Bed_Filtered_VCFs"
ground_truth_snps = "/home/iaymergen/Bed_Filtered_VCFs/filtered_snps.vcf.gz"
ground_truth_indels = "/home/iaymergen/Bed_Filtered_VCFs/filtered_indels.vcf.gz"

# List of files to process
files_to_process = [
    "final_bowtie_mutect_nobase.vcf.recode.vcf",
    "final_bowtie_mutect_Withbase.vcf.recode.vcf",
    "final_bowtie_somaticsniper_no_Base.vcf",
    "final_bowtie_somaticsniper_with_Base.vcf",
    "final_bowtie_strelka_nobase.vcf.recode.vcf",
    "final_bowtie_strelka_Withbase.vcf.recode.vcf",
    "final_bwa_mutect_noBase.vcf.recode.vcf",
    "final_bwa_mutect_WithBase.vcf.recode.vcf",
    "final_bwa_somaticsniper_no_Base.vcf",
    "final_bwa_somaticsniper_withBase.vcf",
    "final_bwa_strelka_noBase.vcf.recode.vcf",
    "final_bwa_strelka_WithBase.vcf.recode.vcf",
]

def load_variants(vcf_path):
    """Load variants from a VCF file into a set."""
    print(f"Loading variants from: {vcf_path}")
    variants = set()
    try:
        with vcfpy.Reader.from_path(vcf_path) as reader:
            for record in reader:
                variants.add((record.CHROM, record.POS, record.REF, tuple(str(alt) for alt in record.ALT)))
    except Exception as e:
        print(f"Error while loading variants from {vcf_path}: {e}")
    print(f"Loaded {len(variants)} variants from {vcf_path}")
    return variants

def calculate_metrics(test_vcf_path, truth_vcf_path):
    """Calculate TP, FP, FN for a test VCF file against the ground truth."""
    print(f"Calculating metrics for: {test_vcf_path}")
    test_variants = load_variants(test_vcf_path)
    truth_variants = load_variants(truth_vcf_path)

    tp = len(test_variants & truth_variants)
    fp = len(test_variants - truth_variants)
    fn = len(truth_variants - test_variants)

    print(f"Metrics for {test_vcf_path}: TP={tp}, FP={fp}, FN={fn}")

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print(f"Precision={precision:.3f}, Recall={recall:.3f}, F1-Score={f1_score:.3f}")
    return tp, fp, fn, precision, recall, f1_score

# Analyze each pipeline for SNPs and Indels
results = []
print("Starting analysis of VCF files...")
for file in files_to_process:
    test_vcf_path = os.path.join(vcf_directory, file)
    print(f"Processing file: {file}")

    # Calculate metrics for SNPs
    print(f"Calculating SNP metrics for {file}")
    snp_tp, snp_fp, snp_fn, snp_precision, snp_recall, snp_f1 = calculate_metrics(test_vcf_path, ground_truth_snps)

    # Calculate metrics for Indels
    print(f"Calculating Indel metrics for {file}")
    indel_tp, indel_fp, indel_fn, indel_precision, indel_recall, indel_f1 = calculate_metrics(test_vcf_path, ground_truth_indels)

    results.append({
        "File": file,
        "SNP_TP": snp_tp, "SNP_FP": snp_fp, "SNP_FN": snp_fn,
        "SNP_Precision": snp_precision, "SNP_Recall": snp_recall, "SNP_F1": snp_f1,
        "Indel_TP": indel_tp, "Indel_FP": indel_fp, "Indel_FN": indel_fn,
        "Indel_Precision": indel_precision, "Indel_Recall": indel_recall, "Indel_F1": indel_f1
    })

# Save results to a CSV
print("Saving results to CSV...")
df = pd.DataFrame(results)
df.to_csv(os.path.join(vcf_directory, "metrics_results_snps_and_indels.csv"), index=False)
print("Metrics saved to metrics_results_snps_and_indels.csv")
