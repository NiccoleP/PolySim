Raw QTL data for height was obtained from the UKBiobank:

wget https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz

The LD blocks from Berisa & Pickrell, 2016:
https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-chr3.bed

The original SLiM demographic model was obtained with the following command. 
python3 -m stdpopsim  -e slim --slim-script --slim-scaling-factor 1 --slim-burn-in 1 -v HomSap -c chr22 -o tree.ts -d OutOfAfrica_3G09 100 100 100 > ORIGINAL.slim

To run the slim script 
slim  -d Q=1.0 -d burn_in=1.0 -d "tree_name='Q1_BN1_tree.ts'" -d "inter_tree_name='pre_Q1_BN1_tree.ts'" -d  mean_p0=-0.03782096 -d mean_p1=-0.0541799 -d mean_p2=-0.05647234 -d sd_p0=0.0112498 -d sd_p1=0.01103384 -d sd_p2=0.01423588 -d "independent_variants='50_irnt_independent_variants'" script.slim & disown
