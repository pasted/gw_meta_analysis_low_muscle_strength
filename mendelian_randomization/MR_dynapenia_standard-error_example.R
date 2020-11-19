## set working directory first due to permission issues
setwd("CHARGE_dynapenia/scripts")

library(data.table)
library(dplyr)
library(tidyverse)
library(readxl)
library(TwoSampleMR)
library(rapportools)
library(Cairo)
library(hash)


outcomes <- c("dynapenia_ewgsop_all","dynapenia_ewgsop_female","dynapenia_ewgsop_male","dynapenia_fnih_all","dynapenia_fnih_female","dynapenia_fnih_male")


sample_size <- c(254894, 134594, 120300, 254894, 134594, 120300)

gwas_path <- "CHARGE_dynapenia/metal_outputs"
run_dir <- "20200313_se_run"

selected_traits <- c("Age lost virginity","Age-related Macular Degeneration","Alcohol consumption","Allergic disease (2018)","Allergic disease (2018)-specific not asthma locus","Alzheimer's Disease (2018 GWAS)","Asthma (2018)","Asthma (2018)-specific not allergic disease locus","Asthma and allergic disease (2018)","Atrial fibrillation","Benign prostatic hyperplasia (2018 GWAS)","Birth weight","BMD (Femoral neck)","BMD (Lumbar spine)","Body Fat Percentage","Breast Cancer (2018 GWAS)","CAD (2015 paper) (no ApoE)","Chronic Kidney Disease (2016 replicated)","Chronotype (23andMe)","Colorectal Cancer (2018 GWAS)","Compound white cell [WBC#] White blood cell count","CRP (2018 GWAS + conditional SNPs)","Dehydroepiandrosterone sulphate (DHEAS)","Depression (2018 meta-analysis)","Educational attainment","eGFR","Fasting glucose","Fasting insulin","Ferritin (log)","FEV1","FEV1/FVC","Fracture risk (2019)","FVC","Glucose challenge (2 hours)","HDL","Homocysteine","Immature red cell [RET#] Reticulocyte count","Insulin resistance","Insulin secretion","Intelligence","Isoleucine","LDL","Liver enzyme (ALP)","Liver enzyme (GGT)","Lymphoid white cell [LYMPH#] Lymphocyte count","Major Depressive Disorder","Mature red cell [HGB] Hemoglobin concentration","Mature red cell [MCHC] Mean corpuscular hemoglobin concentration","Mature red cell [MCV] Mean corpuscular volume","Mature red cell [RBC#] Red blood cell count","Mature red cell [RDW] Red cell distribution width","Menarche (2017 + UKB)","Migraine","Myeloid white cell [BASO#] Basophil count","Myeloid white cell [EO#] Eosinophil count","Myeloid white cell [GRAN#] Granulocyte count","Myeloid white cell [MONO#] Monocyte count","Myeloid white cell [MYELOID#] Myeloid white cell count","Myeloid white cell [NEUT#] Neutrophil count","Osteoarthritis (2018 GWAS)","Osteoarthritis (2019 meta-analysis UKB)","Osteoarthritis (2019 meta-analysis UKB) hip-specific","Osteoarthritis (2019 meta-analysis UKB) knee-specific","Parents survival (2019) GWAS replicated only","Physical activity (MVPA)","Platelet [PCT] Plateletcrit","Platelet [PDW] Platelet distribution width","Platelet [PLT#] Platelet count","Prostate Cancer (2018 paper)","Resting Heart Rate","Rheumatoid Arthritis","Schizophrenia (2018)","Sex hormone-binding globulin (SHBG)","Stroke (any ischemic)","Telomere length","Transferrin","Triglycerides","Type-2 Diabetes (2018 GWAS)","Vitamin D (25(OH)D) concentration (Jiang 2018)","Whole body lean mass","WHR (adjBMI Men)","WHR (adjBMI SexCombined)","WHR (adjBMI Women)")

	
today <- Sys.time()
current_date <- format(today, format="%Y-%m-%d")

results_path <- paste0("CHARGE_dynapenia/results/MR/", current_date)

ifelse(!dir.exists(file.path(results_path)), dir.create(file.path(results_path)), FALSE)

## read in the data trimming white space due to presence in various fields, that are used to match records on


labels <- read.csv2("GRS_labels_v3_20191218.csv", sep=",", strip.white = TRUE, na=c("",NA))

this_data <- read.csv2("SNP_list_191218.edited-headers.txt", sep="\t", strip.white = TRUE, na=c("",NA))


## merge labels

merged_data <- merge(this_data, labels, by.x=c("friendly", "source_units_beta"), by.y=c("friendly", "source_units"), all = TRUE)



## Setup two additional blank columns named as required by TwoSampleMR

merged_data$beta <- NULL
merged_data$se <- NULL

## populate beta.exposure and se.exposure based on whether the external_beta field is present or not
## if present then the exposure beta should be the external beta
## else let the the exposure beta be the UK Biobank beta

merged_data <- as.data.table(merged_data)

merged_data[is.na(external_beta), beta := as.double(as.character(beta_weighting_meta_analysis))]
merged_data[!is.na(external_beta), beta := as.double(as.character(external_beta))]

merged_data[is.na(external_se), se := as.double(as.character(se_for_beta))]
merged_data[!is.na(external_se), se := as.double(as.character(external_se))]


merged_data <- setnames(merged_data, old=c("rsid", "effect_allele_freq"), new=c("SNP", "eaf"), skip_absent=TRUE)

## subset data on the presence of a grs_var code, select only the useful columns

selected_df <- subset(merged_data, (!is.na(grs_var)), select=c(trait, SNP, effect_allele, other_allele, eaf, beta, se, grs_var, do_mr, sub_cats, sex_specific, desc, filename))

## sub_cats unfortunately is a factor, convert to characters to test
## [1] ""                            "Biomarker or risk measure"  
## [3] "Cancer"                      "Hormones"                   
## [5] "Inflammation or Auto-immune" "Lifecourse factors"         
## [7] "Lifestyle risks"             "Lipids"                     
## [9] "Metabolic syndrome"          "Other disease"              
## [11] "Single variant"              "Vitamins"     

## select only traits with do_mr == TRUE - traits with only 1 or 2 SNPs discounted, as are some unpublished GRS
selected_df <- subset(selected_df, do_mr == TRUE)

## remove rows with no effect allele - check the original files
sub_cat_df_clean <- subset(selected_df, !is.na(effect_allele) || effect_allele != "" || effect_allele != "-")

## Only the 83 selected traits
sub_cat_df_clean <- subset(sub_cat_df_clean, trait %in% selected_traits)

## This will be the unique vector of phenotypes to loop over and run the analysis against

traits <- as.character(unique(sub_cat_df_clean$grs_var))

# Loop over traits and subset by that trait

#I:/Projects/Garan/CHARGE_dynapenia/metal_outputs/dynapenia_ewgsop_all/20200215_ss_run/metal_ewgsop_all.ss-se-combined.tsv


for(outcome in outcomes) {
	##remove sex specific GRS if required by outcome cohort
	##if(grepl("female", outcome)  ) {
	##	sub_cat_df_clean <- subset(sub_cat_df_clean, sex_specific != "male")
	##}else if(grepl("male", outcome) ) {
	##	sub_cat_df_clean <- subset(sub_cat_df_clean, sex_specific != "female")
	##}

	outcome_dir <- paste(c(gwas_path, outcome, "20200313_se_run"), collapse="/")
	#metal_ukb_dynapenia-ewgsop-all_2020-03-13.ext-id.gc.hetdf.tsv
	#dynapenia_ewgsop_all_20200313_se_run.ext-id.gc.hetdf.tsv
	
	input_file <- paste0(outcome_dir, "/", outcome, "_", run_dir, ".ext-id.gc.hetdf.tsv")
		
	ifelse(!dir.exists(file.path(results_path, outcome)), dir.create(file.path(results_path, outcome)), FALSE)
	
	setwd(file.path(results_path, outcome))
	
	result_filename <- paste( c(paste(c(outcome, "results", current_date), collapse="_"), "csv"), collapse=".")	
	result_full_path <- file.path(results_path, outcome, result_filename)
	
	reordered_filename <- paste( c(paste(c(outcome, "reordered", current_date), collapse="_"), "csv"), collapse=".")
	
	reordered_full_path <- file.path(results_path, outcome, reordered_filename)
	
	## convert dataframe to exposure
	phenotype_exp_dat <- format_data(sub_cat_df_clean, phenotype_col="trait", type="exposure")
    
	#chromosome      position        ext_id  allele1 allele2 freq1   effect  se      pvalue  rsid
	#1       100000012       1_100000012_G_T T       G       0.2361  -0.0016 9e-4    0.07964 rs10875231
	#1       100000827       1_100000827_C_T T       C       0.2854  -0.0012 9e-4    0.1783  rs6678176

	## Phenotype field used to separate the analysis on an exposure by exposure basis
    outcome_dat <- read_outcome_data(snps = phenotype_exp_dat$SNP,
    filename = input_file,
    sep = "\t",
    snp_col = "snp_id",
    beta_col = "effect",
    se_col = "se",
    effect_allele_col = "allele1",
    other_allele_col = "allele2",
    eaf_col = "freq1",
    pval_col = "pvalue")
    
    outcome_dat$outcome <- outcome
    
	## harmonise data
	
	dat <- harmonise_data(exposure_dat = phenotype_exp_dat, outcome_dat = outcome_dat)
	
	## run MR and save results to dataframe
	
	res <- mr(dat)
	pleiotropy <- mr_pleiotropy_test(dat)
	pleiotropy <- setnames(pleiotropy, old=c("se", "pval"), new=c("pleiotropy_se", "pleiotropy_pval"), skip_absent=TRUE)
	pleiotropy <- subset(pleiotropy, select=c("id.exposure", "egger_intercept","pleiotropy_se", "pleiotropy_pval") )
	
	mr_egger <- subset(res, method=="MR Egger")
	mr_egger <- setnames(mr_egger, old=c("b", "se", "pval"), new=c("mr_egger_b", "mr_egger_se", "mr_egger_pval"), skip_absent=TRUE)
	mr_egger <- subset(mr_egger, select=c("id.exposure", "mr_egger_b", "mr_egger_se", "mr_egger_pval") )
	
	simple <- subset(res, method=="Simple mode")
	simple <- setnames(simple, old=c("b", "se", "pval"), new=c("simple_b", "simple_se", "simple_pval"), skip_absent=TRUE)
	simple <- subset(simple, select=c("id.exposure", "simple_b", "simple_se", "simple_pval") )
	
	weighted <- subset(res, method=="Weighted mode")
	weighted <- setnames(weighted, old=c("b", "se", "pval"), new=c("weighted_b", "weighted_se", "weighted_pval"), skip_absent=TRUE)
	weighted <- subset(weighted, select=c("id.exposure", "weighted_b", "weighted_se", "weighted_pval") )
	
	ivw <- subset(res, method=="Inverse variance weighted")
	ivw <- setnames(ivw, old=c("b", "se", "pval"), new=c("ivw_b", "ivw_se", "ivw_pval"), skip_absent=TRUE)
	ivw <- subset(ivw, select=c("id.exposure", "id.outcome", "outcome", "exposure", "nsnp", "ivw_b", "ivw_se", "ivw_pval") )
	
	reordered_results <- list(ivw, mr_egger, simple, weighted, pleiotropy) %>% reduce(left_join, by="id.exposure")
	
	
	fwrite(res, file = result_full_path)
	fwrite(reordered_results, file = reordered_full_path)
	
	## run MR on all available exposures
	mr_report(dat)
	
}
