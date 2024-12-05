import pandas as pd

input_file = snakemake.input.data
steps_qc = snakemake.params.steps_qc
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)
data = data[data['feature_name']=='steps']

cols = ['subject_id', 'device_id', 'test_run_begin', 'test_repeat', 
        'test_run_type', 'day_of_study', 'numeric_value']
data = data.filter(cols)

rename = {'numeric_value': 'n_steps'}
data = data.rename(columns=rename)
        
data['steps_ok'] = data['n_steps'] <= steps_qc
    
data.to_csv(output_file, index=None, header=True)
