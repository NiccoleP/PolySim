#!/usr/bin/env python
# -*- coding: utf-8 -*-

rule ukbb_sumstats:
    output:
        "ukbiobank/{code}.gwas.imputed_v3.{sex}.tsv.bgz"
    shell:
        "cat ukbiobank/UKBBGWAS_Imputed_v3_File_Manifest_Release_20180731_Manifest_201807.tsv |" \
        " grep -w {wildcards.code} | grep {sex} | awk '{{print $12}}' | xargs wget"

rule ld_blocks:
    output:
        "ld_blocks/fourier_ls-{chr}.bed"
    shell:
        "wget https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-{chr}.bed"
        "echo "$(tail -n +2 fourier_ls-{wildcards.chr}.bed)" > fourier_ls-{wildcards.chr}.bed"

rule variants:
    output:
        "variants/variants.tsv.bgz"
        shell:
        "wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
rule change_R_script:
    output:
        "Poly_Sim_chr{number}.R"
    shell:
       "sed 's/\^3/\^{wilcards.number}/g' Poly_Sim.R | sed 's/{wildcards.chr}/chr{wildcards.number}/g' > Poly_Sim_chr{wildcards.number}"

rule independent_markers:
    input:
        "ukbiobank/{code}.gwas.imputed_v3.{sex}.tsv.bgz",
        "ld_blocks/fourier_ls-chr3.bed"
        "variants/variants.tsv.bgz"
    output:
        "ukbiobank/{code}_independent_variants"
    shell:
        "Rscript Poly_Sim_chr{wildcards.number}.R {input} "
        
rule make_bed: 
    output: 
        "{chr}_variants.bed"
    shell:
        "awk 'BEGIN {{OFS="\t"} {print {wildcards.number},$1-1,$1}}' 50_irnt_independent_variants_chr22 > {wildcards.chr}_variants.bed"
rule demographic_model:
    input:
        "ukbiobank/{code}_independent_variants"
    output:
        "stdpopsim_slim_script_50scaled.slim"
    shell:
         "python3 -m stdpopsim  -e slim --slim-scaling-factor 50 --slim-script --slim-burn-in 7300 -v HomSap -s 1046 -c chr3 -o foo.ts -d OutOfAfrica_3G09 0 1198 0"

