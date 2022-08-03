args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# expr_list.rds
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))
process_kallisto_output(work_dir, 'Kim_kallisto.zip', tx2gene)

