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
        "wget https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-{wildcards.chr}.bed"
        "echo "$(tail -n +2 fourier_ls-{wildcards.chr}.bed)" > fourier_ls-{wildcards.chr}.bed"

rule variants:
    output:
        "variants/variants.tsv.bgz"
        shell:
        "wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"

rule independent_markers:
    input:
        trait="ukbiobank/{code}.gwas.imputed_v3.{sex}.tsv.bgz", 
        blocks="ld_blocks/fourier_ls-{chr}.bed",
        variants="variants/variants.tsv.bgz"
    output:
        "ukbiobank/{code}_independent_variants" 
    shell:
        "Rscript Poly_Sim.R {input.trait} {input.blocks} {input.variants}"
        
rule change_R_script:
    output:
        "Poly_Sim_chr{number}.R"
    shell:
       "sed 's/\^3/\^{wilcards.number}/g' Poly_Sim.R | sed 's/{wildcards.chr}/chr{wildcards.number}/g' > Poly_Sim_chr{wildcards.number}"

rule demographic model:
    output:
        "stdpopsim_slim_script_{scale}_scaled_{burn}_{chr}.slim",
        "tree.ts"
    shell:
        "python3 -m stdpopsim  -e slim --slim-scaling-factor {wildcards.scale} --slim-script --slim-burn-in {wildcards.burn} -v HomSap -c {wildcards.chr} -o {wildcards.tree}.ts -d OutOfAfrica_3G09 {wildcards.YRI} {wildcards.CEU} {wildcards.CHB} \
        > stdpopsim_slim_script_{wildcards.scale}_scaled_{wildcards.burn}_{wildcards.chr}.slim"
