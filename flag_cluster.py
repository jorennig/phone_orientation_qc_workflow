import pandas as pd

input_file = snakemake.input.data[0]
cluster_flag = snakemake.params.cluster_flag
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['cluster_ok'] = ~(data['cluster'].isin(cluster_flag))

data.to_csv(output_file, index=False)
