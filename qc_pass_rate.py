import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

cols = ['subject_id']
data['n_tests'] = data.groupby(cols)['subject_id'].transform('count')

id_vars = ['subject_id', 'device_id', 'test_run_begin', 'test_repeat', 
           'test_run_type', 'day_of_study', 'n_tests']
qc_cols = ['phone_not_on_table_ok', 'span_duration_ok', 'cluster_ok', 
           'steps_ok', 'new_qc_pass']
data = pd.melt(data, id_vars=id_vars, value_vars=qc_cols, 
               var_name='feature_qc', value_name='qc_pass')

group_cols = ['subject_id', 'feature_qc']
data['qc_pass_n'] = data.groupby(group_cols)['qc_pass'].transform('sum')

data['qc_pass_rate_percent'] = data['qc_pass_n'] / data['n_tests'] * 100

drop_cols = ['device_id', 'test_run_begin', 'test_repeat', 
             'test_run_type', 'day_of_study', 'qc_pass', 'qc_pass_n']
data = data.drop(columns=drop_cols) \
           .drop_duplicates()

data['two_week_period'] = 0
data['qc_type'] = 'none'

data.to_csv(output_file, index = None, header = True)
