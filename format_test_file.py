import pandas as pd

input_file = snakemake.input.data
digital_features = snakemake.params.digital_features
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)
data = data[data['feature_name'].isin(digital_features)]

cols = ['subject_id', 'device_id', 'feature_name', 
        'test_run_begin', 'span_begin', 'span_end', 'extract_path',
        'test_repeat', 'test_run_type', 'day_of_study', 'two_week_period', 
        'numeric_value', 'qc_pass']
data = data.filter(cols)

rename = {'feature_name': 'feature_digital', 'numeric_value': 'value_digital'}
data = data.rename(columns=rename)
           
data.to_csv(output_file, index=None, header=True)
