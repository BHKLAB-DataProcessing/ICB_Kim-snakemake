library(data.table)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

## Get Clinical data
clin_original = read.table( file.path(input_dir, "gas_korean_cli_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE , dec=",")
clin_original = clin_original[ clin_original[ , "treatment" ] %in% "pre" , ]
rownames(clin_original) <- clin_original$X

selected_cols <- c("X", "response")
clin = as.data.frame( cbind( clin_original[ , selected_cols ] , "PD-1/PD-L1" , "Gastric" , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA , NA ) )
colnames(clin) = c( "patient" , "recist" , "drug_type" , "primary" , "age" , "histo" , "response" , "pfs" ,"os" , "t.pfs" , "t.os" , "stage" , "sex" , "response.other.info" , "dna" , "rna" )
rownames(clin) = clin$patient

clin$recist = ifelse( clin$recist %in% 0 , "PD" , ifelse( clin$recist %in% -1 , "SD" , ifelse( clin$recist %in% 1 , "CR" , NA ) ) )

clin$response = Get_Response( data=clin )
clin$rna = "fpkm"
clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
expr_samples <- colnames(data.frame(expr_list[['expr_gene_tpm']]))
patient = intersect( expr_samples , rownames(clin) )

for(assay_name in names(expr_list)){
  expr <- data.frame(expr_list[[assay_name]])
  expr =  expr[ , patient ]
  write.table( 
    expr , 
    file= file.path(output_dir, paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')) , 
    quote=FALSE , 
    sep=";" , 
    col.names=TRUE , 
    row.names=TRUE 
  )
}

clin = clin[ patient , ]
clin_original <- clin_original[ patient, ]
clin <- format_clin_data(clin_original, 'X', selected_cols, clin)

# snv = read.csv( file.path(output_dir, "SNV.csv") , sep=";" , stringsAsFactors=FALSE )
# snv = snv[ snv$Sample %in% patient , ]
# snv_patient = unique( sort( snv$Sample ) )

case = cbind( patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = case[,1]
# case[ snv_patient , "snv" ] = 0

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Kim', annotation_tissue=annotation_tissue, check_histo=FALSE)

annotation_drug <- read.csv(file=file.path(annot_dir, 'curation_drug.csv'))
clin <- add_column(clin, treatmentid='', .after='tissueid')

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
