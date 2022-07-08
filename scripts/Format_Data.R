library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")

#############################################################################
#############################################################################

expr = read.table( file.path(input_dir, "gas_korean_exp_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE )
rownames(expr) = expr[ , 1 ]
expr = expr[ , -1 ]

#############################################################################
#############################################################################
## Get Clinical data

clin_original = read.table( file.path(input_dir, "gas_korean_clin_data.csv") , sep="," , header=TRUE , stringsAsFactors = FALSE , dec=",")
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


#############################################################################
#############################################################################

patient = intersect( colnames(expr) , rownames(clin) )
clin = clin[ patient , ]
clin_original <- clin_original[ patient, ]
clin <- format_clin_data(clin_original, 'X', selected_cols, clin)

expr =  expr[ , patient ]

snv = read.csv( file.path(output_dir, "SNV.csv") , sep=";" , stringsAsFactors=FALSE )
snv = snv[ snv$Sample %in% patient , ]
snv_patient = unique( sort( snv$Sample ) )

case = cbind( patient , 0 , 0 , 1 )
colnames(case ) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = case[,1]
case[ snv_patient , "snv" ] = 1

write.table( case , file = file.path(output_dir, "cased_sequenced.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( clin , file = file.path(output_dir, "CLIN.csv") , sep = ";" , quote = FALSE , row.names = FALSE)
write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
