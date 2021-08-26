""" evip (Expression-based Varient Impact Phenotyping) process handler """

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

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


def calculate_wildtype_distribution(df: pd.DataFrame, rep_null_dist:list, wt_name:str='RNF43_WT')->dict:
    """ eVIP buildWT_dict """

    df_wt = df.loc[[x.startswith(wt_name) for x in df.index],[x.startswith(wt_name) for x in df.columns]]
    np.fill_diagonal(df_wt.values, np.nan)

    rep_rankpts = df_wt.apply(lambda x: np.percentile(x.dropna().values, 50), axis=1).values
    wt_rep_rankpt = np.percentile(rep_rankpts, 50)

    wt_rep_pval = wilcoxon(rep_rankpts, [rep_null_dist[0] for _ in rep_rankpts]).pvalue

    ret1 = {wt_name: {'wt_rep': wt_rep_rankpt, 'wt_rep_dist': rep_rankpts, 'wt_rep_pval':wt_rep_pval}}
    ret2 = [wt_rep_pval]
    ret3 = [wt_name]

    return ret1, ret2, ret3
