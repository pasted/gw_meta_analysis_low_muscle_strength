#Additonal filtering and QC
#Example for EWGSOP combined dataset

check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c("ggplot2", "readr", "tidyverse", "gdata", "dplyr", "tidyr", "data.table", "stringr", "scales", "CMplot")
check.packages(packages)


dir <- "CHARGE_dynapenia/QC"


out_dir <- "CHARGE_dynapenia/results/run_2020-03-13_1600"
description <- "metal_ukb_dynapenia-ewgsop-all_2020-03-13"

setwd(out_dir)

file_name_metal <- paste0(out_dir, "/METAL-CHARGE_dynapenia_EWGSOP-ALL-2020-03-13-GC.se.1.tbl")

file_name_ukb_snps <- "CHARGE_dynapenia/raw_data/ukbb_ewgsop_all.tsv"

metal <- read_delim(file_name_metal, "\t", escape_double = FALSE, col_types = cols(HetDf = col_integer()), trim_ws = TRUE)

ukb <- read_delim(file_name_ukb_snps, "\t", escape_double = FALSE, col_types = cols(chr = col_integer(), position = col_integer()), trim_ws = TRUE)

metal <- separate(metal, MarkerName, c("chromosome","position"), sep="_", remove=FALSE, extra="drop")


metal <- setnames(metal, "MarkerName", "ext_id")
ukb <- setnames(ukb, "id", "ext_id")

drops <- c("af_coded_all", "oevar_imp", "beta","se","p_value")

ukb <- ukb[ , !(names(ukb) %in% drops)]

#convert both types of dataframe to data.table for compatability
metal.dt <- data.table(metal)
ukb.dt <- data.table(ukb)

metal_with_ukb_rsid <- merge(x=metal, y=ukb, by.x="ext_id", by.y="ext_id", all=FALSE)

metal_with_ukb_rsid  <- setnames(metal_with_ukb_rsid, c("position.x"), c("position"), skip_absent=TRUE )


#keep: 
#Chromosome
#Position
#rsid
#Pvalue
#Allele1
#Allele2
#Effect#
#StdErr
#remove: "alleleA", "alleleB", "MinFreq", "MaxFreq", "FreqSE"
#
#drops <- c("alleleA", "alleleB", "MinFreq","MaxFreq","FreqSE")
#
#metal_with_ukb_rsid_cleaned <- metal_with_ukb_rsid[ , !(names(metal_with_ukb_rsid) %in% drops)]


#Variants with no genomic position: remove them
#1802479	1	NA	1_NA_A_G	a	g	0.9953	0	0.6871	1	NA
#3782774	2	NA	2_NA_A_G	a	g	0.0967	0	0.1353	1	NA
#10079359	6	NA	6_NA_C_T	t	c	0.9852	0	0.3304	1	NA
#10079360	6	NA	6_NA_A_G	a	g	0.0017	0	0.8602	1	NA
#11412281	7	NA	7_NA_A_C	a	c	0.0029	0	0.8627	1	NA
#13723706	9	NA	9_NA_A_T	a	t	0.8956	0	0.1315	1	NA
#13723707	9	NA	9_NA_C_T	t	c	0.0047	0	0.6824	1	NA
#14881633	10	NA	10_NA_A_T	a	t	0.981	0	0.3468	1	NA


metal_ukb_rsid <- metal_with_ukb_rsid %>% drop_na(position)

#if position column has been imported as characters check the string "NA"
metal_ukb_rsid <- subset(metal_ukb_rsid, position != "NA")

#uppercase the METAL alleles

metal_ukb_rsid <- data.table(metal_ukb_rsid)

metal_ukb_rsid$Allele1 <- toupper(metal_ukb_rsid$Allele1)
metal_ukb_rsid$Allele2 <- toupper(metal_ukb_rsid$Allele2)


out_file_name <- paste(description, "ext-id.gc.full.tsv", sep=".")
out_file <- file.path(out_dir, out_file_name)


write_tsv(metal_ukb_rsid, out_file)



metal_ukb_rsid_maf <- subset(metal_ukb_rsid, Freq1 <= 0.99 & Freq1 >= 0.01)

#Output subset of SNPs present in at least 2 cohorts; HetDF = number of cohorts - 1
metal_ukb_hetdf <- subset(metal_ukb_rsid_maf, HetDf > 0)

#nrow(metal_ukb_hetdf)
#[1] 2925497

#remove dashes from "p-value", lowercase & snakecase all column names

metal_ukb_rsid_maf <- setnames(metal_ukb_rsid_maf, c("Allele1","Allele2","Freq1","Zscore","P-value"), c("allele1","allele2","freq1","z_score","pvalue"))

#metal_ukb_rsid_maf <- setnames(metal_ukb_rsid_maf, old=c("marker_name"), new=c("ext_id"))
	
#duplicate rows removed
#nrow = 10123338 at start
metal_ukb_rsid_maf <- distinct(metal_ukb_rsid_maf)
#nrow = 10069649

#rows containing duplicate ext_id remove - keeping first occurrence
metal_ukb_rsid_maf <- metal_ukb_rsid_maf %>% distinct(ext_id, .keep_all=TRUE)
#nrow = 10023922


out_file_name <- paste(description, "ext-id.gc.0-01.tsv", sep=".")
out_file <- file.path(out_dir, out_file_name)


write_tsv(metal_ukb_rsid_maf, out_file)



drops <- c("HetPVal", "HetChiSq", "HetISq", "HetDf", "Direction")
dropcols <- names(metal_ukb_hetdf)[names(metal_ukb_hetdf) %in% drops]
metal_ukb_hetdf[, c(dropcols) := NULL]

#remove dashes from "p-value", lowercase & snakecase all column names
 
metal_ukb_hetdf <- setnames(metal_ukb_hetdf, c("Allele1","Allele2","Freq1","Effect","StdErr","P-value", "alleleA", "alleleB"), c("allele1","allele2","freq1","effect","se","pvalue","effect_allele","noneffect_allele"), skip_absent=TRUE)


metal_ukb_hetdf <- distinct(metal_ukb_hetdf)
#nrow(metal_ukb_hetdf)
#[1] 2905671

#rows containing duplicate ext_id remove - keeping first occurrence
metal_ukb_hetdf <- metal_ukb_hetdf %>% distinct(ext_id, .keep_all=TRUE)

# nrow(metal_ukb_hetdf)
#[1] 2859778



out_file_name <- paste(description, "ext-id.gc.hetdf.tsv", sep=".")
out_file <- file.path(out_dir, out_file_name)


write_tsv(metal_ukb_hetdf, out_file)
