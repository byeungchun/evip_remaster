""" kallisto handler """

import numpy as np
import pandas as pd

def concat_tpm_from_abundance(abundance_files:dict)->pd.DataFrame:
    """ read tpm from abundence files and concatenate all """

    tpm_df = list()
    for sample, _file in abundance_files.items():
        df = pd.read_csv(_file, sep='\t', index_col=0)
        df = df.tpm
        df.name = sample
        tpm_df.append(df)
    
    return pd.DataFrame(tpm_df).T

def sum_tpm_by_gene_name(df_tpm:pd.DataFrame, df_gtf:pd.DataFrame)->pd.DataFrame:
    """ summerise tpm values by gene name """

    df_gtf = df_gtf[['transcript_id','gene_name']]
    df_gtf = df_gtf[df_gtf.transcript_id.str.strip() != ''].drop_duplicates()

    df_tpm.loc[:,'transcript_id'] = [x.split('.')[0] for x in df_tpm.index]

    df = pd.merge(df_tpm, df_gtf, on='transcript_id')
    df = df[df.columns.delete(-2)] # remove transcript_id column
    _columns = [x for x in set(df.gene_name) if x != '']

    df = df.groupby('gene_name').agg(sum).loc[_columns,:]

    return df

def calculate_tpm_log2(df:pd.DataFrame, min_log_val:float=1.0)->pd.DataFrame:
    """ calculate tpm log2 values """

    return df[df >= min_log_val].apply(np.log2).dropna(how='all')

def calculate_tpm_zscore(df:pd.DataFrame)-> pd.DataFrame:
    """ calculate z-score values of tpm log2 dataframe """

    return df.apply(lambda x: (x - df.mean(axis=1)) / df.std(axis=1)).fillna(0.0)

def calculate_corr_sample(df:pd.DataFrame)->pd.DataFrame:
    """ calculate correlation of RNA samples """

    return df.corr()