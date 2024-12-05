import pandas as pd
from scipy.stats import spearmanr

def correlation(chunk, value_cols):
    n = chunk.shape[0]
    r, p = spearmanr(chunk[value_cols[0]], 
                     chunk[value_cols[1]], 
                     nan_policy='omit')
    return pd.Series([r, p, n], index=['rho', 'p', 'n'])

input_file = snakemake.input.data[0]
group_cols = snakemake.params.group_cols
value_cols = snakemake.params.value_cols
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})

corr = data.groupby(group_cols) \
           .apply(correlation, value_cols) \
           .reset_index()

corr.to_csv(output_file, index=None, header=True)
