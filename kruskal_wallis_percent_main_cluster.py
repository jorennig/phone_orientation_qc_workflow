import pandas as pd
from scipy.stats import kruskal

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'main_cluster': str}
data = pd.read_csv(input_file, dtype=dtype)

labels = sorted(list(data['main_cluster'].unique()))
value_array = []
for l in labels:
    cluster = data.loc[data['main_cluster']==l, 'changes_percent_cluster']
    value_array.append(cluster)
s, p = kruskal(*value_array)

kw_res = pd.DataFrame([s, p], index=['stat', 'p'])

kw_res.to_csv(output_file, index=True, header=True)
