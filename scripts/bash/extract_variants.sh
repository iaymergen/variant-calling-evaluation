#!/bin/bash

# Create a unified variant list
output_variant_list="all_variants.txt"
> $output_variant_list

# Extract variants from each VCF file
for vcf in *.vcf; do
    echo "Processing $vcf"
    grep -v '^#' $vcf | awk '{print $1"\t"$2"\t"$4"\t"$5}' >> $output_variant_list
done

# Remove duplicates
sort $output_variant_list | uniq > $output_variant_list

echo "Unique variants saved in $output_variant_list"

