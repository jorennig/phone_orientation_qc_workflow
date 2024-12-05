import pandas as pd

input_file = snakemake.input.data[0]
feature = snakemake.params.feature
value = snakemake.params.value
cluster_col = snakemake.params.cluster_col
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

summary = data.groupby([feature, cluster_col])[value] \
              .agg(n='count', median='median', mad='mad') \
              .reset_index()

summary.to_csv(output_file, index = None, header = True)
