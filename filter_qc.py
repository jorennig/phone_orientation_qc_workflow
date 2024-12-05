import pandas as pd

input_file = snakemake.input.data[0]
qc_variable = snakemake.params.qc_variable
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data = data[data[qc_variable]].drop(columns=[qc_variable])

data.to_csv(output_file, index=False)
