import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_rel, ttest_ind


# Load the CSV file
metrics_file = "/home/iaymergen/Bed_Filtered_VCFs/metrics_results_snps_and_indels.csv"
metrics_data = pd.read_csv(metrics_file)

# Add additional columns for configurations
metrics_data["Mapper"] = metrics_data["File"].apply(lambda x: "bowtie" if "bowtie" in x else "bwa")
metrics_data["Caller"] = metrics_data["File"].apply(
    lambda x: "mutect" if "mutect" in x else "somaticsniper" if "somaticsniper" in x else "strelka"
)
metrics_data["Recalibration"] = metrics_data["File"].apply(lambda x: "WithBase" if "WithBase" in x else "NoBase")

# Calculate Accuracy for SNPs and Indels
metrics_data["SNP_Accuracy"] = metrics_data.apply(
    lambda row: row["SNP_TP"] / (row["SNP_TP"] + row["SNP_FP"] + row["SNP_FN"]) if (row["SNP_TP"] + row["SNP_FP"] + row["SNP_FN"]) > 0 else 0,
    axis=1
)
metrics_data["Indel_Accuracy"] = metrics_data.apply(
    lambda row: row["Indel_TP"] / (row["Indel_TP"] + row["Indel_FP"] + row["Indel_FN"]) if (row["Indel_TP"] + row["Indel_FP"] + row["Indel_FN"]) > 0 else 0,
    axis=1
)

# 1. Heatmaps for SNP and Indel F1-Scores
for variant_type in ["SNP", "Indel"]:
    heatmap_data = metrics_data.pivot_table(values=f"{variant_type}_F1", index="Mapper", columns="Caller")
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".4f")
    plt.title(f"{variant_type} F1-Scores by Mapper and Caller")
    plt.tight_layout()
    plt.show()

# 2. Heatmaps for Accuracy (SNP and Indel)
for variant_type in ["SNP", "Indel"]:
    heatmap_data = metrics_data.pivot_table(values=f"{variant_type}_Accuracy", index="Mapper", columns="Caller")
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".4f")
    plt.title(f"{variant_type} Accuracy by Mapper and Caller")
    plt.tight_layout()
    plt.show()

# 3. Histograms for Metrics (Precision, Recall, F1, Accuracy)
metrics = ["Precision", "Recall", "F1", "Accuracy"]
for metric in metrics:
    for variant_type in ["SNP", "Indel"]:
        plt.figure(figsize=(10, 6))
        sns.histplot(data=metrics_data, x=f"{variant_type}_{metric}", hue="Caller", kde=True, palette="muted", bins=20)
        plt.title(f"Distribution of {variant_type} {metric} Across Callers")
        plt.xlabel(f"{variant_type} {metric}")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.show()

# 4. Boxplots for Metrics (Precision, Recall, F1, Accuracy)
for metric in metrics:
    for variant_type in ["SNP", "Indel"]:
        plt.figure(figsize=(10, 6))
        sns.boxplot(x="Caller", y=f"{variant_type}_{metric}", hue="Recalibration", data=metrics_data, palette="viridis")
        plt.title(f"{variant_type} {metric} by Caller and Recalibration")
        plt.ylabel(f"{variant_type} {metric}")
        plt.xlabel("Variant Caller")
        plt.tight_layout()
        plt.show()

# 5. Venn Diagram (replace with real data if available)
plt.figure(figsize=(8, 8))
venn3(subsets=(100, 80, 60, 30, 40, 20, 10), set_labels=("Bowtie", "BWA", "Strelka"))
plt.title("Overlap of Variants Called by Pipelines")
plt.tight_layout()
plt.show()

# 6. PCA Plot
features = ["SNP_F1", "Indel_F1", "SNP_Precision", "SNP_Recall", "SNP_Accuracy", "Indel_Accuracy"]
pca_data = metrics_data[features].dropna()

# Standardize data
scaler = StandardScaler()
pca_data_scaled = scaler.fit_transform(pca_data)

# Apply PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pca_data_scaled)
metrics_data["PCA1"] = pca_result[:, 0]
metrics_data["PCA2"] = pca_result[:, 1]

plt.figure(figsize=(10, 8))
sns.scatterplot(x="PCA1", y="PCA2", hue="Caller", style="Recalibration", data=metrics_data, palette="deep")
plt.title("PCA Plot of Pipeline Metrics")
plt.tight_layout()
plt.show()

# 7. Statistical Testing for SNP and Indel Metrics
for metric in metrics:
    for variant_type in ["SNP", "Indel"]:
        with_base = metrics_data[metrics_data["Recalibration"] == "WithBase"][f"{variant_type}_{metric}"].reset_index(drop=True)
        no_base = metrics_data[metrics_data["Recalibration"] == "NoBase"][f"{variant_type}_{metric}"].reset_index(drop=True)

        print(f"\nStatistical Testing for {variant_type} {metric}:")
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
