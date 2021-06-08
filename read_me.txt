Raw QTL data for height was obtained from the UKBiobank:

wget https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz

The LD blocks from Berisa & Pickrell, 2016:
https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-chr3.bed

The original SLiM demographic model was obtained with the following command. 
python3 -m stdpopsim  -e slim --slim-script --slim-scaling-factor 1 --slim-burn-in 1 -v HomSap -c chr22 -o tree.ts -d OutOfAfrica_3G09 100 100 100 > ORIGINAL.slim
