import pandas as pd

input_file = snakemake.input.data[0]
span_duration_bounds = snakemake.params.span_duration_bounds
output_file = snakemake.output[0]

date_columns = ['span_begin', 'span_end']
dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype, parse_dates=date_columns)

data['span_duration'] = (data['span_end'] - data['span_begin']).dt.total_seconds()
data['span_duration_ok'] = data['span_duration'].between(*span_duration_bounds)

data.to_csv(output_file, index=False)
