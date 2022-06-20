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
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Kim-data/main/"
patients = pd.read_csv(data_source + "annot_WES.txt", sep="\t", header=0)["patient"].values

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    shell:
        """
        Rscript -e \
        '
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Kim", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule format_data:
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv")
    input:
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "download/gas_korean_clin_data.csv"),
        S3.remote(prefix + "download/gas_korean_exp_data.csv")
    shell:
        """
        Rscript scripts/Format_Data.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_snv:
    output:
        S3.remote(prefix + "processed/SNV.csv")
    input: 
        S3.remote(prefix + "download/annot_WES.txt"),
        S3.remote(expand(prefix + "download/annot_vcf/{patient}.txt", patient=patients.sort()))
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_annot_vcf:
    output:
        S3.remote(expand(prefix + "download/annot_vcf/{patient}.txt", patient=patients.sort()))
    input:
        S3.remote(prefix + "download/annot_WES.txt")
    run:
        for p in patients:
            shell("wget {data_source}annot_vcf/{p}.txt -O {prefix}download/annot_vcf/{p}.txt")

rule download_data:
    output:
        S3.remote(prefix + "download/annot_WES.txt"),
        S3.remote(prefix + "download/gas_korean_clin_data.csv"),
        S3.remote(prefix + "download/gas_korean_exp_data.csv")
    shell:
        """
        wget {data_source}annot_WES.txt -O {prefix}download/annot_WES.txt
        wget {data_source}gas_korean_clin_data.csv -O {prefix}download/gas_korean_clin_data.csv
        wget {data_source}gas_korean_exp_data.csv -O {prefix}download/gas_korean_exp_data.csv
        """ 