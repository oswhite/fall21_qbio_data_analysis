# QBIO Mid Semester Project
# Olivia White
# Do the most commonly mutated genes in gastric cancer have any correlation with young or old patients?


library (BiocManager)
library (TCGAbiolinks)
library (maftools)

#load data
maf_file <- data.table::fread("C:/Users/livis/Documents/fall21_qbio_data_analysis/week5_maf/GDCdata/TCGA.STAD.mutect.c06465a3-50e7-46f7-b2dd-7bd654ca206b.DR-10.0.somatic.maf.csv") #Replace this with appropriate file path and name

clinic <- read.csv("/Users/livis/Documents/fall21_qbio_data_analysis/coad_clinical_data.csv", row.names = 1)
colnames( clinic )[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_dataframe <- read.maf(maf_file, isTCGA = TRUE, clinicalData = clinic)

#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)


#generate MAF summary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


#generate oncoplot of data
oncoplot( maf_dataframe, top=5 )


#generate subsets for age
bool_young <-  clinic$age_at_initial_pathologic_diagnosis < 45 
young_patients <- clinic$Tumor_Sample_Barcode[ bool_young ]
young_maf <- subsetMaf( maf_dataframe, tsb = young_patients )

bool_old <- clinic$age_at_initial_pathologic_diagnosis >= 45
old_patients <- clinic$Tumor_Sample_Barcode[bool_old]
old_maf <- subsetMaf( maf_dataframe, tsb = old_patients)


num_young <- length (young_patients)
num_old <- length (old_patients)


#generate subset for TTN mutations
TTN_maf <- subsetMaf( maf_dataframe, genes = "TTN" )

#generate numbers for contingency table
TTN_mutated_samples <- TTN_maf@clinical.data$Tumor_Sample_Barcode
num_TTN_mutated_samples <- length( TTN_mutated_samples )

TTN_and_young_samples <- intersect( TTN_mutated_samples, young_patients )
num_TTN_and_young_samples <- length( TTN_and_young_samples )

num_TTN_mutation_old <- num_TTN_mutated_samples - num_TTN_and_young_samples
num_young_No_TTN_mutation <- num_young - num_TTN_and_young_samples

total_samples <- length( maf_dataframe@clinical.data$Tumor_Sample_Barcode )
num_no_TTN_old <- total_samples - num_TTN_and_young_samples - num_TTN_mutation_old - num_young_No_TTN_mutation

#generate contingency table for ttn
contig_table_TTN <- matrix(c(num_TTN_and_young_samples, num_young_No_TTN_mutation,num_TTN_mutation_old, num_no_TTN_old), nrow=2)
contig_table_TTN

#generate Fisher's Exact results for ttn

fe_results_TTN <- fisher.test(contig_table_TTN)
fe_results_TTN


#generate subset for tp53 mutations
TP53_maf <- subsetMaf( maf_dataframe, genes = "TP53" )

#generate numbers for contingency table
TP53_mutated_samples <- TP53_maf@clinical.data$Tumor_Sample_Barcode
num_TP53_mutated_samples <- length( TP53_mutated_samples )

TP53_and_young_samples <- intersect( TP53_mutated_samples, young_patients )
num_TP53_and_young_samples <- length( TP53_and_young_samples )

num_TP53_mutation_old <- num_TP53_mutated_samples - num_TP53_and_young_samples
num_young_No_TP53_mutation <- num_young - num_TP53_and_young_samples

total_samples <- length( maf_dataframe@clinical.data$Tumor_Sample_Barcode )
num_no_TP53_old <- total_samples - num_TP53_and_young_samples - num_TP53_mutation_old - num_young_No_TP53_mutation

#generate contingency table for tp53
contig_table_tp53 <- matrix(c(num_TP53_and_young_samples, num_young_No_TP53_mutation,num_TP53_mutation_old, num_no_TP53_old), nrow=2)
contig_table_tp53

#generate Fisher's Exact results for tp53

fe_results_tp53 <- fisher.test(contig_table_tp53)
fe_results_tp53


#generate subset for MUC16 mutations
MUC16_maf <- subsetMaf( maf_dataframe, genes = "MUC16" )

#generate numbers for contingency table
MUC16_mutated_samples <- MUC16_maf@clinical.data$Tumor_Sample_Barcode
num_MUC16_mutated_samples <- length( MUC16_mutated_samples )

MUC16_and_young_samples <- intersect( MUC16_mutated_samples, young_patients )
num_MUC16_and_young_samples <- length( MUC16_and_young_samples )

num_MUC16_mutation_old <- num_MUC16_mutated_samples - num_MUC16_and_young_samples
num_young_No_MUC16_mutation <- num_young - num_MUC16_and_young_samples

total_samples <- length( maf_dataframe@clinical.data$Tumor_Sample_Barcode )
num_no_MUC16_old <- total_samples - num_MUC16_and_young_samples - num_MUC16_mutation_old - num_young_No_MUC16_mutation

#generate contingency table for MUC16
contig_table_MUC16 <- matrix(c(num_MUC16_and_young_samples, num_young_No_MUC16_mutation,num_MUC16_mutation_old, num_no_MUC16_old), nrow=2)
contig_table_MUC16

#generate Fisher's Exact results for MUC16

fe_results_MUC16 <- fisher.test(contig_table_MUC16)
fe_results_MUC16


#generate subset for LRP1B mutations
LRP1B_maf <- subsetMaf( maf_dataframe, genes = "LRP1B" )

#generate numbers for contingency table
LRP1B_mutated_samples <- LRP1B_maf@clinical.data$Tumor_Sample_Barcode
num_LRP1B_mutated_samples <- length( LRP1B_mutated_samples )

LRP1B_and_young_samples <- intersect( LRP1B_mutated_samples, young_patients )
num_LRP1B_and_young_samples <- length( LRP1B_and_young_samples )

num_LRP1B_mutation_old <- num_LRP1B_mutated_samples - num_LRP1B_and_young_samples
num_young_No_LRP1B_mutation <- num_young - num_LRP1B_and_young_samples

total_samples <- length( maf_dataframe@clinical.data$Tumor_Sample_Barcode )
num_no_LRP1B_old <- total_samples - num_LRP1B_and_young_samples - num_LRP1B_mutation_old - num_young_No_LRP1B_mutation

#generate contingency table for LRP1B
contig_table_LRP1B <- matrix(c(num_LRP1B_and_young_samples, num_young_No_LRP1B_mutation,num_LRP1B_mutation_old, num_no_LRP1B_old), nrow=2)
contig_table_LRP1B

#generate Fisher's Exact results for LRP1B

fe_results_LRP1B <- fisher.test(contig_table_LRP1B)
fe_results_LRP1B


#generate subset for SYNE1 mutations
SYNE_maf <- subsetMaf( maf_dataframe, genes = "SYNE1" )

#generate numbers for contingency table
SYNE_mutated_samples <- SYNE_maf@clinical.data$Tumor_Sample_Barcode
num_SYNE_mutated_samples <- length( SYNE_mutated_samples )

SYNE_and_young_samples <- intersect( SYNE_mutated_samples, young_patients )
num_SYNE_and_young_samples <- length( SYNE_and_young_samples )

num_SYNE_mutation_old <- num_SYNE_mutated_samples - num_SYNE_and_young_samples
num_young_No_SYNE_mutation <- num_young - num_SYNE_and_young_samples

total_samples <- length( maf_dataframe@clinical.data$Tumor_Sample_Barcode )
num_no_SYNE_old <- total_samples - num_SYNE_and_young_samples - num_SYNE_mutation_old - num_young_No_SYNE_mutation

#generate contingency table for syne
contig_table_syne <- matrix(c(num_SYNE_and_young_samples, num_young_No_SYNE_mutation,num_SYNE_mutation_old, num_no_SYNE_old), nrow=2)
contig_table_syne

#generate Fisher's Exact results for syne

fe_results_syne <- fisher.test(contig_table_syne)
fe_results_syne



