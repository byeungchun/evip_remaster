""" evip (Expression-based Varient Impact Phenotyping) process handler """

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, kruskal
from statsmodels.stats.multitest import multipletests

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
    wt_rep_rankpt, rep_rankpts = _calculate_self_connectivity(df_wt)

    wt_rep_pval = wilcoxon(rep_rankpts, [rep_null_dist[0] for _ in rep_rankpts]).pvalue

    ret1 = {wt_name: {'wt_rep': wt_rep_rankpt, 'wt_rep_dist': rep_rankpts, 'wt_rep_pval':wt_rep_pval}}
    ret2 = [wt_rep_pval]
    ret3 = [wt_name]

    return ret1, ret2, ret3

def calculate_wt_mt_comparison(df: pd.DataFrame, rep_null_dist:list, conn_null_dist:list, wt_dict: dict, sig_id:dict)->dict:
    """ eVIP Build comparison part """

    # mut_rep_pvals = list()
    # mut_wt_conn_pvals = list()
    # mut_wt_rep_pvals = list()
    # mut_wt_rep_vs_wt_mut_conn_pvals = list()

    wt_name = [k for k,v in sig_id.items() if v['type']=='WT'][0]
    ret_val = {wt_name: dict()}

    for cell_key, cell_val in sig_id.items():
        if cell_val['type'] != 'MT':
            continue
        
        #STEP1: calculate self connectivity and p-value
        df_mt = df.loc[[x.startswith(cell_key) for x in df.index],[x.startswith(cell_key) for x in df.columns]]
        mut_rankpt, mut_rankpt_dist = _calculate_self_connectivity(df_mt)

        # wilcoxon test on Scipy does not generate a same result as in R if sample is too small
        self_pval = wilcoxon(mut_rankpt_dist, [rep_null_dist[0] for _ in mut_rankpt_dist]).pvalue
        #mut_rep_pvals.append(self_pval)

        #STEP2: calculate connectivity and p-value
        df_mt_wt = df.loc[[x.startswith(wt_name) for x in df.index],[x.startswith(cell_key) for x in df.columns]]
        mut_wt_conn_dist = list(df_mt_wt.apply(lambda x: np.percentile(x.values, 50)).values)
        mut_wt_conn_dist.extend(list(df_mt_wt.apply(lambda x: np.percentile(x.values, 50), axis=1).values))
        mut_wt_conn_rankpt = np.percentile(mut_wt_conn_dist, 50)
        conn_pval = wilcoxon(mut_wt_conn_dist, [conn_null_dist[0] for _ in mut_wt_conn_dist]).pvalue
        #mut_wt_conn_pvals.append(conn_pval)

        #STEP3: connectivity between mutant self connectivity and wildtype dist.
        mut_wt_rep_pval = wilcoxon(mut_rankpt_dist, wt_dict[wt_name]['wt_rep_dist']).pvalue
        #mut_wt_rep_pvals.append(mut_wt_rep_pval)

        #STEP4: Kruskal test
        wt_mut_rep_vs_wt_mut_conn_pval = kruskal(wt_dict[wt_name]['wt_rep_dist'], mut_rankpt_dist, mut_wt_conn_dist).pvalue
        #mut_wt_rep_vs_wt_mut_conn_pvals.append(wt_mut_rep_vs_wt_mut_conn_pval)

        _median = [np.median(mut_rankpt_dist), np.median(mut_wt_conn_rankpt), np.median(wt_dict[wt_name]['wt_rep_dist'])]
        median_diff = max(_median) - min(_median)

        ret_val[wt_name][cell_key] = {
            'wt_rep': rep_null_dist[0],
            'mut_rep': mut_rankpt,
            'mut_rep_pval': self_pval,
            'mut_wt_connectivity': mut_wt_conn_rankpt,
            'mut_wt_conn_pval': conn_pval,
            'mut_wt_rep_pval': mut_wt_rep_pval,
            'wt_mut_rep_vs_wt_mut_conn_pval':wt_mut_rep_vs_wt_mut_conn_pval,
            'kruskal_diff': median_diff
        }

    mut_wt_rep_pvals = [mut_val['mut_wt_rep_pval'] for mut_val in ret_val[wt_name].values()]
    mut_wt_conn_pvals = [mut_val['mut_wt_conn_pval'] for mut_val in ret_val[wt_name].values()]
    mut_wt_rep_vs_wt_mut_conn_pvals = [mut_val['wt_mut_rep_vs_wt_mut_conn_pval'] for mut_val in ret_val[wt_name].values()]

    #STEP5: calculate corrected pvalues
    for i, mut_val in enumerate(ret_val[wt_name].values()):
        mut_val.update({'mut_wt_rep_c_pvals': multipletests(mut_wt_rep_pvals, method='fdr_bh')[1][i]})
        mut_val.update({'mut_wt_conn_c_pvals': multipletests(mut_wt_conn_pvals, method='fdr_bh')[1][i]})
        mut_val.update({'mut_wt_rep_vs_wt_mut_conn_c_pvals': multipletests(mut_wt_rep_vs_wt_mut_conn_pvals, method='fdr_bh')[1][i]})

    return ret_val

def predict_mutant_functionality(
    res_comp:dict,
    use_c_pval:bool=True,
    mut_wt_thresh:float=0.1,
    mut_wt_rep_diff:float=0.0,
    c_thresh:float=0.1,
    disting_thresh:float=0.1,
    cond_median_max_diff_thresh:float=0.2
)->dict:
    """ eVIP_predict function """

    wt_name = list(res_comp.keys())[0]
    
    # for mut_key, mut_val in res_comp[wt_name].items():
    for mut_key, mut_val in res_comp[wt_name].items():
        wt_rep = mut_val['wt_rep']
        mut_rep = mut_val['mut_rep']
        mut_wt_conn = mut_val['mut_wt_connectivity']
        kruskal_diff = mut_val['kruskal_diff']
        if use_c_pval:
            mut_wt_rep_pval = mut_val['mut_wt_rep_c_pvals']
            mut_wt_conn_pval = mut_val['mut_wt_conn_c_pvals']
            disting_pval = mut_val['wt_mut_rep_vs_wt_mut_conn_c_pvals']
        else:
            mut_wt_rep_pval = mut_val['mut_wt_rep_pval']
            mut_wt_conn_pval = mut_val['mut_wt_conn_pval']
            disting_pval = mut_val['wt_mut_rep_vs_wt_mut_conn_pval']


        res_comp[wt_name][mut_key]['prediction'] = _get_prediction_6(
            wt_rep,
            mut_rep,
            mut_wt_conn,
            mut_wt_rep_pval,
            mut_wt_conn_pval,
            disting_pval,
            mut_wt_thresh,
            mut_wt_rep_diff,
            c_thresh,
            disting_thresh,
            cond_median_max_diff_thresh,
            kruskal_diff
        )

    return res_comp
    
def _get_prediction_6(
    wt_rep,
    mut_rep,
    mut_wt_conn,
    mut_wt_rep_pval,
    mut_wt_conn_pval,
    disting_pval,
    mut_wt_thresh,
    mut_wt_rep_diff,
    c_thresh,
    disting_thresh,
    cond_median_max_diff_thresh,
    cond_median_max_diff
)->str:
    """ Determin prediction """
    
    ret_val = 'NI'
    if disting_pval < disting_thresh:
        if _max_diff(wt_rep, mut_rep, mut_wt_conn) < mut_wt_rep_diff:
            ret_val = 'Neutral' if mut_wt_conn_pval < c_thresh else 'Error'

        if cond_median_max_diff > cond_median_max_diff_thresh:
            ret_val = 'COF'
            if mut_wt_rep_pval < mut_wt_thresh:
                if wt_rep < mut_rep:
                    ret_val = 'GOF' if mut_rep - wt_rep >= mut_wt_rep_diff else 'COF'
                elif wt_rep > mut_rep:
                    ret_val = 'LOF' if wt_rep - mut_rep >= mut_wt_rep_diff else 'COF'

    if mut_wt_conn_pval < c_thresh:
        ret_val = 'Neutral'

    return ret_val

def _max_diff(wt_rep, mut_rep, mut_wt_conn):
    max_diff = abs(wt_rep - mut_rep)

    wt_conn_diff = abs(wt_rep - mut_wt_conn)
    if wt_conn_diff > max_diff:
        max_diff = wt_conn_diff

    mut_conn_diff = abs(mut_rep - mut_wt_conn)
    if mut_conn_diff > max_diff:
        max_diff = mut_conn_diff

    return max_diff
    

def _calculate_self_connectivity(df:pd.DataFrame)->list:
    """ get median values """
    
    np.fill_diagonal(df.values, np.nan)
    rep_rankpts = df.apply(lambda x: np.percentile(x.dropna().values, 50), axis=1).values
    wt_rep_rankpt = np.percentile(rep_rankpts, 50)

    return wt_rep_rankpt, rep_rankpts