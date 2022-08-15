import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]

rule get_MultiAssayExp:
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv"),
        # S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb=3000,
        disk_mb=3000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v40.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Kim", input_dir = paste0("{prefix}", "processed"), expr_with_counts_isoforms=TRUE), 
            "{prefix}{filename}"
        );
        '
        """

rule format_data:
    input:
        # S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "download/gas_korean_cli_data.csv"),
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv")
    shell:
        """
        Rscript scripts/Format_Data.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation
        """

# rule format_snv:
#     output:
#         S3.remote(prefix + "processed/SNV.csv")
#     input: 
#         S3.remote(prefix + "download/annot_WES.txt"),
#         S3.remote(prefix + "download/annot_vcf.zip")
#     resources:
#         mem_mb=3000
#     shell:
#         """
#         unzip -d {prefix}download/ {prefix}/download/annot_vcf.zip
#         Rscript scripts/Format_SNV.R \
#         {prefix}download \
#         {prefix}processed \
#         """

rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/Kim_kallisto.zip"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + "download/expr_list.rds")
    shell:
        """
        Rscript scripts/format_downloaded_data.R \
        {prefix}download \
        {prefix}annotation 
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv 
        """

rule download_data:
    output:
        S3.remote(prefix + "download/Kim_kallisto.zip"),
        S3.remote(prefix + "download/gas_korean_cli_data.csv")
    shell:
        """
        wget -O {prefix}download/gas_korean_cli_data.csv https://github.com/xmuyulab/ims_gene_signature/raw/main/data/gas_korean_cli_data.csv
        wget -O {prefix}download/Kim_kallisto.zip https://zenodo.org/record/6968411/files/Kim_kallisto.zip?download=1
        """ 