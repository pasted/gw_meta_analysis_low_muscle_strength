#!/bin/bash

## annotate variables with paths to your version of UKB + phenotypes
## syntax is illustrative to show the exact options/parameters used

## run BOLT-LMM on UK Biobank sample
BOLT-LMM_v2.3.2/bolt \
    --bed=${BFILE_IMP_BED} \
    --bim=${BFILE_IMP_BIM} \
    --fam=${BFILE_IMP_FAM} \
    --phenoFile=${FILE} \
    --phenoCol=${PHENO} \
    --modelSnps=${MODEL_SNPS} \
    --lmm \
    --LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --statsFile=${OUT_FILE_GENO} \
    --bgenFile=${BGEN} \
    --sampleFile=${SAMPLE} \
    --noBgenIDcheck \
    --statsFileBgenSnps=${OUT_FILE_IMPUTED} \
    --numThreads 8 \
    --bgenMinMAF=0.0002 \
    --verboseStats \
    --lmmForceNonInf \
    --covarFile=${COVAR_FILE} \
    --covarMaxLevels 30 \
    --qCovarCol=age --covarCol=sex --covarCol=axiom_bileve --covarCol=assessment_centre

