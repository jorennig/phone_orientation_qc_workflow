import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

summary = data['cluster'].value_counts() \
                         .reset_index(name='n') \
                         .rename(columns={'index': 'cluster'})

summary.to_csv(output_file, index=None, header=True)
