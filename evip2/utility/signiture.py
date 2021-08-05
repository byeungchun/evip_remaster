""" Signiture file handling """

import pandas as pd

from itertools import chain

def extract_sample_from_sigfile(
    sigfile:str,
    delim:str = r'\t',
    sample_delim:str='|',
    sample_colname:str='distil_id')->list:
    """ extract sample id from signiture file """

    df = pd.read_csv(sigfile, sep=delim)

    return list(chain.from_iterable(df[sample_colname].apply(lambda x: x.split(sample_delim))))
