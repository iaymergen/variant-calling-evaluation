import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel, ttest_ind

# 1. Load the CSV file
metrics_file = "/home/iaymergen/Bed_Filtered_VCFs/metrics_results_snps_and_indels.csv"
metrics_data = pd.read_csv(metrics_file)

# 2. Add additional columns for configurations
metrics_data["Mapper"] = metrics_data["File"].apply(lambda x: "bowtie" if "bowtie" in x else "bwa")
metrics_data["Caller"] = metrics_data["File"].apply(
    lambda x: "mutect" if "mutect" in x else "somaticsniper" if "somaticsniper" in x else "strelka"
)
metrics_data["Recalibration"] = metrics_data["File"].apply(lambda x: "WithBase" if "WithBase" in x else "NoBase")

# 3. Calculate Metrics for SNPs and Indels
for variant_type in ["SNP", "Indel"]:
    # Calculate Precision
    metrics_data[f"{variant_type}_Precision"] = metrics_data.apply(
        lambda row: row[f"{variant_type}_TP"] / (row[f"{variant_type}_TP"] + row[f"{variant_type}_FP"]) if (row[f"{variant_type}_TP"] + row[f"{variant_type}_FP"]) > 0 else 0,
        axis=1
    )
    # Calculate Recall
    metrics_data[f"{variant_type}_Recall"] = metrics_data.apply(
        lambda row: row[f"{variant_type}_TP"] / (row[f"{variant_type}_TP"] + row[f"{variant_type}_FN"]) if (row[f"{variant_type}_TP"] + row[f"{variant_type}_FN"]) > 0 else 0,
        axis=1
    )
    # Calculate F1-Score
    metrics_data[f"{variant_type}_F1"] = metrics_data.apply(
        lambda row: 2 * (row[f"{variant_type}_Precision"] * row[f"{variant_type}_Recall"]) / (row[f"{variant_type}_Precision"] + row[f"{variant_type}_Recall"]) if (row[f"{variant_type}_Precision"] + row[f"{variant_type}_Recall"]) > 0 else 0,
        axis=1
    )
    # Calculate Accuracy
    metrics_data[f"{variant_type}_Accuracy"] = metrics_data.apply(
        lambda row: row[f"{variant_type}_TP"] / (row[f"{variant_type}_TP"] + row[f"{variant_type}_FP"] + row[f"{variant_type}_FN"]) if (row[f"{variant_type}_TP"] + row[f"{variant_type}_FP"] + row[f"{variant_type}_FN"]) > 0 else 0,
        axis=1
    )

# 4. Summarize Metrics
grouped = metrics_data.groupby(["Caller", "Recalibration"]).mean(numeric_only=True)
print("Mean Metrics by Caller and Recalibration:")
print(grouped[[f"{variant_type}_{metric}" for variant_type in ["SNP", "Indel"] for metric in ["Precision", "Recall", "F1", "Accuracy"]]])

# 5. Visualizations

## a. Compare F1-Scores for SNPs and Indels
melted_data = metrics_data.melt(
    id_vars=["File", "Mapper", "Caller", "Recalibration"],
    value_vars=["SNP_F1", "Indel_F1"],
    var_name="Metric_Type",
    value_name="F1-Score"
)

plt.figure(figsize=(14, 8))
sns.barplot(x="Caller", y="F1-Score", hue="Recalibration", data=melted_data, palette="viridis")
plt.title("SNP and Indel F1-Scores Across Variant Callers and Recalibration")
plt.ylabel("F1-Score")
plt.xlabel("Variant Caller")
plt.legend(title="Recalibration")
plt.tight_layout()
plt.show()

## b. Heatmaps for Metrics
for metric in ["Precision", "Recall", "F1", "Accuracy"]:
    heatmap_data = metrics_data.pivot_table(values=f"SNP_{metric}", index="Caller", columns="Recalibration")
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="YlGnBu", fmt=".4f")
    plt.title(f"SNP {metric} by Variant Caller and Recalibration")
    plt.tight_layout()
    plt.show()

## c. Boxplots for Metrics
for metric in ["Precision", "Recall", "F1", "Accuracy"]:
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Caller", y=f"SNP_{metric}", hue="Recalibration", data=metrics_data, palette="coolwarm")
    plt.title(f"SNP {metric} by Caller and Recalibration")
    plt.ylabel(f"SNP {metric}")
    plt.xlabel("Variant Caller")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Caller", y=f"Indel_{metric}", hue="Recalibration", data=metrics_data, palette="coolwarm")
    plt.title(f"Indel {metric} by Caller and Recalibration")
    plt.ylabel(f"Indel {metric}")
    plt.xlabel("Variant Caller")
    plt.tight_layout()
    plt.show()

# 6. Statistical Testing for SNPs
for metric in ["F1", "Precision", "Recall", "Accuracy"]:
    with_base = metrics_data[metrics_data["Recalibration"] == "WithBase"][f"SNP_{metric}"].reset_index(drop=True)
    no_base = metrics_data[metrics_data["Recalibration"] == "NoBase"][f"SNP_{metric}"].reset_index(drop=True)

    print(f"\nStatistical Testing for SNP {metric}:")
    if len(with_base) > 1 and len(no_base) > 1:
        if len(with_base) == len(no_base):
            # Paired t-test
            t_stat, p_value = ttest_rel(with_base, no_base)
            print(f"Paired T-test: t_stat={t_stat:.4f}, p_value={p_value:.4f}")
        else:
            # Independent t-test
            min_length = min(len(with_base), len(no_base))
            t_stat, p_value = ttest_ind(with_base[:min_length], no_base[:min_length], equal_var=False)
            print(f"Independent T-test: t_stat={t_stat:.4f}, p_value={p_value:.4f}")
    else:
        print("Not enough data points for statistical testing.")
