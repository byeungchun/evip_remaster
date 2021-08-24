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

def calculate_null_distribution(df: pd.DataFrame, sig_id:dict)->dict:
    """ eVIP getNullDist function """

    rand_gfp = 'GFP_4'
    base_gfp = 'GFP_1'

    #STEP1: Control cell - MT, WT distribution
    conn_null_dist = list()
    for _cell in sig_id.keys() - ['GFP']:
        conn_null_dist.append(np.percentile(df.loc[rand_gfp, df.columns.str.contains(_cell)], 50))
    
    #STEP2: Control cell distribution
    rep_null_dist = list()
    gfp_cols = [x for x in df.columns if x.startswith('GFP') and x != base_gfp]
    for _ in range(len(sig_id.keys() - ['GFP'])):
        gfp_vals = df.loc[df.index.str.contains(base_gfp), gfp_cols].stack().values
        rep_null_dist.append(np.percentile(gfp_vals,50))
    
    return rep_null_dist, conn_null_dist


