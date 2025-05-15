#!/bin/bash

# Directory to store extracted variant files
output_dir="variants_output"
mkdir -p $output_dir

# List of files to process
files=(
    "final_bowtie_mutect_no_Base_onlyPass.vcf"
    "final_bowtie_mutect_withBase_onlyPass.vcf"
    "final_bowtie_somaticsniper_no_Base.vcf"
    "final_bowtie_somaticsniper_with_Base.vcf"
    "final_bowtie_strelka_no_Base_onlyPass.vcf"
    "final_bowtie_strelka_with_Base_onlyPass.vcf"
    "final_bwa_mutect_noBase_onlyPass.vcf"
    "final_bwa_mutect_withBase_onlyPass.vcf"
    "final_bwa_somaticsniper_no_Base.vcf"
    "final_bwa_somaticsniper_with_Base.vcf"
    "final_bwa_strelka_no_Base_onlyPass.vcf"
    "final_bwa_strelka_with_Base_onlyPass.vcf"
    "final_galaxy_bowtie_strelka_onlyPass.vcf"
    "final_galaxy_bwa_strelka_onlyPass.vcf"
)

# Loop through each file and extract variants
for file in "${files[@]}"; do
    if [[ -f $file ]]; then
        echo "Processing $file"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $file > $output_dir/variants_${file%%.*}.txt
    elif [[ -f $file.gz ]]; then
        echo "Processing compressed file $file.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $file.gz > $output_dir/variants_${file%%.*}.txt
    else
        echo "File $file or $file.gz not found, skipping."
    fi
done

echo "Variant extraction completed. Results saved in $output_dir."

