library(data.table)
library(readxl) 
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_abundance.R')

# EXPR.txt
dir.create(file.path(work_dir, 'rnaseq'))
unzip(file.path(work_dir, 'Kim_kallisto.zip'), exdir=file.path(work_dir, 'rnaseq'))

samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names=FALSE, recursive = FALSE)

expr <- process_abundance(file.path(work_dir, 'rnaseq'), samples)

write.table(expr, file.path(work_dir, 'EXPR.tsv'), col.names = TRUE, row.names = TRUE, sep = '\t')
