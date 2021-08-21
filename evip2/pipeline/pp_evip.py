""" pipeline default eVIP """

from gtfparse import read_gtf

from evip2.utility.signiture import extract_sample_from_sigfile
from evip2.utility.kallisto import (
    concat_tpm_from_abundance,
    sum_tpm_by_gene_name,
    calculate_tpm_log2,
    calculate_tpm_zscore
)
from evip2.utility.evip import calculate_corr_sample

def main(
    sigfile:str,
    gtffile:str,
    abundance_folder:str
    )->bool:
    """ eVIP runner """

    # STEP1: read signiture files
    samples = extract_sample_from_sigfile(sigfile)

    # STEP2: load tpm from kallisto abundence files
    abundance_files = {x:f'{abundence_folder}/{x}/abundance.tsv' for x in samples}
    df_tpm = concat_tpm_from_abundance(abundance_files)

    # STEP3: sum tpm values by gene names
    df_gtf = read_gtf(gtffile)
    df_gene_tpm = sum_tpm_by_gene_name(df_tpm, df_gtf)

    # STEP4: calculate log2 of gene tpm table
    df_gene_tpm_log = calculate_tpm_log2(df_gene_tpm)

    # STEP5: calculate z-score of gene tpm log table
    df_gene_tpm_zscore = calculate_tpm_zscore(df_gene_tpm_log)

    # STEP6: calculate RNA sample correlation
    df_sample_corr = calculate_corr_sample(df_gene_tpm_zscore)
