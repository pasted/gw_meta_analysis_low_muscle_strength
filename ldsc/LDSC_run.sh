#!/bin/bash

## load python environment on cluster
module load Anaconda2/4.2.0
source activate ldsc

## use LDSC munge function to make sure summary stats are formatted correctly
IFS='/' read -ra FILES <<< "dynapenia_ewgsop_all.Jones_2020.for_ldsc.txt/Gripmax.UKB_2019.for_ldsc.txt/Height.Yengo_2018.Wood_plus_UKB.for_ldsc.txt/BMI.Yengo_2018.Locke_plus_UKB.for_ldsc.txt/WHRadjBMI.Pulit_2019.for_ldsc.txt/Lean_mass_appendicular.Zillikens_2017.for_ldsc.txt/Lean_mass_wholebody.Zillikens_2017.for_ldsc.txt/AD_sumstats_Jansenetal.txt/CAD.van_der_Haarst_2018.for_ldsc.txt/CKD_Pattaro_2016.txt/Mahajan.NatGenet2018b.T2D.European.txt/Michailidou_BC_2017.txt/Osteoarthritis_Zengini_2018.txt/RA.Okada_2014.for_ldsc.txt/Schumacher_PC_2018.txt/Stroke_EUR_Malik_2018.txt/Fracture.Trajanoska_2018.for_ldsc.txt"
N_FILES=${#FILES[@]}
for (( i=0; i<${N_FILES}; i++ ));
do

	ldsc/munge_sumstats.py \
		--sumstats ${FILES[$i]} \
		--out ${FILES[$i]} \
		--merge-allelesldsc/eur_w_ld_chr/w_hm3.snplist

done

## run LDSC on list of GWAS to get genetic correlations with low grip strength (dynapenia)
ldsc/ldsc.py \
	--rg data/dynapenia_ewgsop_all.Jones_2020.for_ldsc.txt.sumstats.gz,data/Gripmax.UKB_2019.for_ldsc.txt.sumstats.gz,data/Height.Yengo_2018.Wood_plus_UKB.for_ldsc.txt.sumstats.gz,data/BMI.Yengo_2018.Locke_plus_UKB.for_ldsc.txt.sumstats.gz,data/WHRadjBMI.Pulit_2019.for_ldsc.txt.sumstats.gz,data/Lean_mass_appendicular.Zillikens_2017.for_ldsc.txt.sumstats.gz,data/Lean_mass_wholebody.Zillikens_2017.for_ldsc.txt.sumstats.gz,data/AD_sumstats_Jansenetal.txt.sumstats.gz,data/CAD.van_der_Haarst_2018.for_ldsc.txt.sumstats.gz,data/CKD_Pattaro_2016.txt.sumstats.gz,data/Mahajan.NatGenet2018b.T2D.European.txt.sumstats.gz,data/Michailidou_BC_2017.txt.sumstats.gz,data/Osteoarthritis_Zengini_2018.txt.sumstats.gz,data/RA.Okada_2014.for_ldsc.txt.sumstats.gz,data/Schumacher_PC_2018.txt.sumstats.gz,data/Stroke_EUR_Malik_2018.txt.sumstats.gz,data/Fracture.Trajanoska_2018.for_ldsc.txt.sumstats.gz,data/parents_mortality.Timmers_2019.txt.sumstats.gz \
	--ref-ld-chr ldsc/eur_w_ld_chr/ \
	--w-ld-chr 7ldsc/eur_w_ld_chr/ \
	--out lowgrip_LDSC_genetic_correlations


