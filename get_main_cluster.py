import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

cols = ['subject_id']
data['main_cluster'] = data.groupby(cols)['cluster'] \
                           .transform(lambda x: x.value_counts().idxmax())

data['main_cluster_n'] = data.groupby(cols)['cluster'] \
                             .transform(lambda x: x.value_counts().max())

data['percent_main_cluster'] = data['main_cluster_n'] / data['n_tests'] * 100

data.to_csv(output_file, index=None, header=True)
