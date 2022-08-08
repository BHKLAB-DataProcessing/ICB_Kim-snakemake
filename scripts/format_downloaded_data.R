library(tximport)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# expr_list.rds
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))

dir.create(file.path(work_dir, 'rnaseq'))
unzip(file.path(work_dir, 'Kim_kallisto.zip'), exdir=file.path(work_dir, 'rnaseq'))
unlink(file.path(work_dir, 'rnaseq', '__MACOSX'), recursive = TRUE)

process_kallisto_output(work_dir, tx2gene)

