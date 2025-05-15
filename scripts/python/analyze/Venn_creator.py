from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Load the variants from the extracted files (example paths)
with open("mutect_variants.txt") as f:
    mutect_variants = set(line.strip() for line in f)
with open("strelka_variants.txt") as f:
    strelka_variants = set(line.strip() for line in f)
with open("somaticsniper_variants.txt") as f:
    somaticsniper_variants = set(line.strip() for line in f)

# Calculate overlaps
only_mutect = mutect_variants - strelka_variants - somaticsniper_variants
only_strelka = strelka_variants - mutect_variants - somaticsniper_variants
only_somaticsniper = somaticsniper_variants - mutect_variants - strelka_variants
mutect_strelka = mutect_variants & strelka_variants - somaticsniper_variants
mutect_somaticsniper = mutect_variants & somaticsniper_variants - strelka_variants
strelka_somaticsniper = strelka_variants & somaticsniper_variants - mutect_variants
all_three = mutect_variants & strelka_variants & somaticsniper_variants

# Create Venn diagram
venn3(subsets=(len(only_mutect), len(only_strelka), len(mutect_strelka),
               len(only_somaticsniper), len(mutect_somaticsniper), len(strelka_somaticsniper),
               len(all_three)),
      set_labels=('Mutect', 'Strelka', 'SomaticSniper'))

# Customize and show plot
plt.title("Overlap of Variants Called by Mutect, Strelka, and SomaticSniper")
plt.show()
