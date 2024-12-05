import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

#data[['x_org', 'y_org', 'z_org']] = data[['x', 'y', 'z']]
data[['x', 'y', 'z']] = abs(data[['x', 'y', 'z']])

data.to_csv(output_file, index=None, header=True)
