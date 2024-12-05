import pandas as pd
from scipy.stats import kruskal

def kruskal_wallis(chunk):
    labels = sorted(list(chunk[cluster_col].unique()))
    value_array = []
    for l in labels:
        cluster = chunk.loc[chunk[cluster_col]==l, value]
        value_array.append(cluster)
    s, p = kruskal(*value_array)
    return pd.Series([s, p], index=['stat', 'p'])

input_file = snakemake.input.data[0]
feature = snakemake.params.feature
value = snakemake.params.value
cluster_col = snakemake.params.cluster_col
output_file = snakemake.output[0]

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)

kw_res = data.groupby(feature).apply(kruskal_wallis).reset_index()

kw_res.to_csv(output_file, index=False, header=True)
