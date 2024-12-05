import pandas as pd
from scipy.stats import mannwhitneyu

def u_test(chunk):
    x = chunk.loc[chunk['qc_value']=='True', 'value_digital']
    y = chunk.loc[chunk['qc_value']=='False', 'value_digital']
    u, p = mannwhitneyu(x, y)
    return pd.Series([u, p, len(x), len(y)], index=['u', 'p', 'n_pass', 'n_flagged'])

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'qc_value': str}
data = pd.read_csv(input_file, dtype=dtype)

cols = ['qc_type', 'feature_digital']
u_test = data.groupby(cols).apply(u_test).reset_index()

u_test.to_csv(output_file, index = None, header = True)
