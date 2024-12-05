import pandas as pd
import numpy as np

def find_max_ratio(row):
    v = sorted([row['x'], row['y'], row['z']], reverse=True)
    max_ratio = v[0]/v[1]
    return max_ratio

input_file = snakemake.input.data[0]
max_ratio = snakemake.params.max_ratio
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['max_ratio'] = data.apply(find_max_ratio, axis=1)
data['definite_orientation'] = data['max_ratio'] >= max_ratio
data['max_ratio_log'] = np.log10(data['max_ratio'])

data.to_csv(output_file, index=None, header=True)
