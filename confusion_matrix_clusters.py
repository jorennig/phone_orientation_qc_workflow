import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str, 'cluster': str})

table = pd.crosstab(data['main_orientation'], data['cluster'])

table.to_csv(output_file, index=True, header=True)
