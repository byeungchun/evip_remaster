""" evip test """

import pandas as pd
from evip2.utility.evip import (
    calculate_null_distribution,
    calculate_wildtype_distribution,
    calculate_wt_mt_comparison
)

sig_id = {
    'GFP':{'type':'CT'},
    'RNF43_G659fs':{'type':'MT'},
    'RNF43_R117fs':{'type':'MT'},
    'RNF43_WT':{'type':'WT'}
}

def test_calculate_null_distribution():
    """ test evip calculate_null_distribution """

    df = pd.read_json(r'./test/data/cell_corr_tbl.json')

    rep_null_dist, conn_null_dist = calculate_null_distribution(df, sig_id)

    assert len(conn_null_dist) == len(rep_null_dist)
    assert rep_null_dist[0] == rep_null_dist[1]
    assert conn_null_dist[0] != conn_null_dist[1]

def test_calculate_wildtype_distribution():
    """ test wildtype distribution """

    df = pd.read_json(r'./test/data/cell_corr_tbl.json')
    rep_null_dist, conn_null_dist = calculate_null_distribution(df, sig_id)

    ret1, ret2, ret3 = calculate_wildtype_distribution(df, rep_null_dist, 'RNF43_WT')

    assert isinstance(ret1, dict)
    assert type(ret2) == type(ret3)
    assert len(ret1['RNF43_WT'].keys()) == 3

def test_calculate_wt_mt_comparison():
    """ test calculate wt mt comparison """

    df = pd.read_json(r'./test/data/cell_corr_tbl.json')
    rep_null_dist, conn_null_dist = calculate_null_distribution(df, sig_id)
    wt_dict, wt_rep_pvals, wt_allele_ordered = calculate_wildtype_distribution(df, rep_null_dist, 'RNF43_WT')

    ret_comp = calculate_wt_mt_comparison(df, rep_null_dist, conn_null_dist, wt_dict, sig_id)

    assert True