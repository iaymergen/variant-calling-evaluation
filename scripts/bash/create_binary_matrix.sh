#!/bin/bash

# Input files
variant_list="all_variants.txt"
output_matrix="binary_matrix.tsv"

# Add header row with VCF file names
echo -e "Variant\t$(ls *.vcf | sed 's/.vcf//g' | tr '\n' '\t')" > $output_matrix

# Process each variant
while read -r chrom pos ref alt; do
    variant="${chrom}_${pos}_${ref}_${alt}"
    binary_row="$variant"
    for vcf in *.vcf; do
        if bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf | grep -q -P "^${chrom}\t${pos}\t${ref}\t${alt}$"; then
            binary_row="${binary_row}\t1"
        else
            binary_row="${binary_row}\t0"
        fi
    done
    echo -e "$binary_row" >> $output_matrix
done < $variant_list

echo "Binary matrix saved to $output_matrix"

