import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel, ttest_ind

# Load the TXT file
metrics_file = "/home/ege/Desktop/Bed_Filtered_VCFs/Pass/gz/metrics_summary1.txt"
metrics_data = pd.read_csv(metrics_file, sep="\t")

# Add additional columns for configurations
metrics_data["Mapper"] = metrics_data["FILE"].apply(lambda x: "bowtie" if "bowtie" in x else "bwa")
metrics_data["Caller"] = metrics_data["FILE"].apply(
    lambda x: "mutect" if "mutect" in x else "somaticsniper" if "somaticsniper" in x else "strelka"
)
metrics_data["Recalibration"] = metrics_data["FILE"].apply(
    lambda x: "WithBase" if "with_Base" in x or "withBase" in x else "NoBase"
)

# 1. Heatmaps for SNP F1-Scores
heatmap_data = metrics_data.pivot_table(values="F1", index="Mapper", columns="Caller")
plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".4f")
plt.title("SNP F1-Scores by Mapper and Caller")
plt.tight_layout()
plt.show()

# 2. Heatmaps for Accuracy (SNP only)
heatmap_data = metrics_data.pivot_table(values="ACCURACY", index="Mapper", columns="Caller")
plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".4f")
plt.title("SNP Accuracy by Mapper and Caller")
plt.tight_layout()
plt.show()

# 3. Histograms for Metrics (Precision, Recall, F1, Accuracy)
metrics = ["PRECISION", "RECALL", "F1", "ACCURACY"]
for metric in metrics:
    plt.figure(figsize=(10, 6))
    sns.histplot(data=metrics_data, x=metric, hue="Caller", kde=True, palette="muted", bins=20)
    plt.title(f"Distribution of SNP {metric} Across Callers")
    plt.xlabel(f"SNP {metric}")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()

# 4. Boxplots for Metrics (Precision, Recall, F1, Accuracy)
for metric in metrics:
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="Caller", y=metric, hue="Recalibration", data=metrics_data, palette="viridis")
    plt.title(f"SNP {metric} by Caller and Recalibration")
    plt.ylabel(f"SNP {metric}")
    plt.xlabel("Variant Caller")
    plt.tight_layout()
    plt.show()

# 5. PCA Plot
features = ["F1", "PRECISION", "RECALL", "ACCURACY"]
pca_data = metrics_data[features].dropna()

# Standardize data
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
pca_data_scaled = scaler.fit_transform(pca_data)

# Apply PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(pca_data_scaled)
metrics_data["PCA1"] = pca_result[:, 0]
metrics_data["PCA2"] = pca_result[:, 1]

plt.figure(figsize=(10, 8))
sns.scatterplot(x="PCA1", y="PCA2", hue="Caller", style="Recalibration", data=metrics_data, palette="deep")
plt.title("PCA Plot of Pipeline Metrics")
plt.tight_layout()
plt.show()

# 6. Statistical Testing for SNP Metrics
for metric in metrics:
    with_base = metrics_data[metrics_data["Recalibration"] == "WithBase"][metric].reset_index(drop=True)
    no_base = metrics_data[metrics_data["Recalibration"] == "NoBase"][metric].reset_index(drop=True)

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
