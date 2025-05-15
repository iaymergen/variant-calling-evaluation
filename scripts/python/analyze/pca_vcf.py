import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load the binary matrix
binary_matrix = pd.read_csv('binary_matrix.tsv', sep='\t', index_col=0)

# Fill NaN values with 0
binary_matrix = binary_matrix.fillna(0)

# Perform PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(binary_matrix.T)  # Transpose to make pipelines as rows
explained_variance = pca.explained_variance_ratio_

# Extract PCA weights (loadings) for variants
variant_weights = pd.DataFrame(pca.components_.T, 
                                index=binary_matrix.index,  # Variants
                                columns=['PC1', 'PC2'])

# Save variant weights to a file
variant_weights.to_csv('variant_pca_weights.tsv', sep='\t')

# Plot PCA results for pipelines
plt.figure(figsize=(10, 8))
for i, pipeline in enumerate(binary_matrix.columns):
    plt.scatter(pca_result[i, 0], pca_result[i, 1], label=pipeline)

plt.title('PCA of VCF Pipelines')
plt.xlabel(f'PC1 ({explained_variance[0]*100:.2f}% variance)')
plt.ylabel(f'PC2 ({explained_variance[1]*100:.2f}% variance)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid()
plt.tight_layout()
plt.savefig('pca_vcf_visualization.png', dpi=300)
plt.show()

# Plot PCA weights for variants in PC1
plt.figure(figsize=(10, 6))
variant_weights['PC1'].sort_values(ascending=False).head(20).plot(kind='bar', color='blue')
plt.title('Top 20 Variant Contributions to PC1')
plt.ylabel('Weight')
plt.xlabel('Variants')
plt.tight_layout()
plt.savefig('variant_weights_PC1.png', dpi=300)
plt.show()

# Plot PCA weights for variants in PC2
plt.figure(figsize=(10, 6))
variant_weights['PC2'].sort_values(ascending=False).head(20).plot(kind='bar', color='green')
plt.title('Top 20 Variant Contributions to PC2')
plt.ylabel('Weight')
plt.xlabel('Variants')
plt.tight_layout()
plt.savefig('variant_weights_PC2.png', dpi=300)
plt.show()

