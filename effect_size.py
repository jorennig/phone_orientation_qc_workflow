import pandas as pd
from numpy import std, mean, sqrt
from itertools import combinations

def cohens_d(chunk, cluster_col, value_col):    
    cluster_labels = chunk[cluster_col].unique()
    
    x_selection = chunk[cluster_col]==cluster_labels[0]
    y_selection = ~(x_selection)
    x = chunk[x_selection][value_col]
    y = chunk[y_selection][value_col]
    
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    d = abs(mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)
    return d

input_file = snakemake.input.data[0]
cluster = snakemake.params.cluster
feature = snakemake.params.feature
value = snakemake.params.value
output_file = snakemake.output[0]

dtype = {'subject_id': str, cluster: str}
data = pd.read_csv(input_file, dtype=dtype)

clusters = list(data[cluster].unique())
cluster_pairs = list(combinations(clusters, 2))

effect_size_list = []
for p in cluster_pairs:
        
    data_c = data[data[cluster].isin(p)]
    
    effect_size = data_c.groupby(feature) \
                        .apply(cohens_d, cluster, value) \
                        .reset_index(name='d')
    effect_size['group1'] = p[0]
    effect_size['group2'] = p[1]
    
    effect_size_list.append(effect_size)

effect_size_tot = pd.concat(effect_size_list) \
                    .sort_values(by=[feature, 'group1'])

effect_size_tot.to_csv(output_file, index=None, header=True)
