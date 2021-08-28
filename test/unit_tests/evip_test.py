""" evip test """

import pandas as pd
from evip2.utility.evip import calculate_null_distribution

def test_calculate_null_distribution():
    """ test evip calculate_null_distribution """

    sig_id = {
        'GFP':{'type':'CT'},
        'RNF43_G659fs':{'type':'MT'},
        'RNF43_R117fs':{'type':'MT'},
        'RNF43_WT':{'type':'WT'}
    }

    df = pd.read_json(r'./test/data/cell_corr_tbl.json')

    rep_null_dist, conn_null_dist = calculate_null_distribution(df, sig_id)

    assert len(conn_null_dist) == len(rep_null_dist)
    assert rep_null_dist[0] == rep_null_dist[1]
    assert conn_null_dist[0] != conn_null_dist[1]