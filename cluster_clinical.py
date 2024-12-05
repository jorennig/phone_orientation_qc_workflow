import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data = data[data['feature_orientation']=='main_cluster'] \
           .rename(columns={'value_orientation': 'main_cluster'}) \
           .drop(columns='feature_orientation')

data.to_csv(output_file, index = None, header = True)
