import pandas as pd
import numpy as np

input_file = snakemake.input.data[0]
quantile = snakemake.params.quantile
output_file = snakemake.output[0]

dtype = {'subject_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['n_tests_cutoff'] = np.quantile(data['n_tests'], quantile)
data['n_tests_ok'] = data['n_tests'] >= data['n_tests_cutoff']

data.to_csv(output_file, index = None, header = True)
