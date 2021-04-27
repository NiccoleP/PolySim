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

rule independent_markers:
    input:
        "ukbiobank/{code}.gwas.imputed_v3.{sex}.tsv.bgz",
        "ld_blocks/fourier_ls-chr3.bed"
        "variants/variants.tsv.bgz"
    output:
        "ukbiobank/{code}_independent_variants"
        "summary_stats"
    shell:
        "Rscript Poly_Sim.R {input} "
