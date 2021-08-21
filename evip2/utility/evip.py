""" evip (Expression-based Varient Impact Phenotyping) process handler """

import numpy as np
import pandas as pd

def calculate_corr_sample(df:pd.DataFrame)->pd.DataFrame:
    """ calculate correlation of RNA samples 
    
    Args:
        df: z-score value of gene by cell 
    
    Returns:
        DataFrame: correlation of cell by cell
    """

    return df.corr()

def compare_functional_impact(df:pd.DataFrame, sig_id:dict)->pd.DataFrame:
    """ compare WT and mutant functional impact using cell correlation 

    Args:
        df: correlation table of cells
        sig_id: signature id dict - [cell type]:[sig id]

    Returns:
        DataFrame:
    """

    for key, val in sig_id.items():
        if val == 'MT':
            df_mt = df.loc[[True if x.find(key)>=0 else False for x in df.index],:].T
            df_mt = df_mt.loc[[True if x.find(key)>=0 else False for x in df_mt.index],:]
            np.fill_diagonal(df_mt.values, None)
            pct50_all = np.percentile(df_mt.stack().values, 50)
            pct50_exp = list(df_mt.apply(lambda x: np.percentile(x.dropna(), 50)))
            
            #TODO: getPairwiseComparisons
