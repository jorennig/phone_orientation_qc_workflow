import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['steps_ok'] = data['steps_ok'].astype(int)

data['two_week_period'] = data['day_of_study'] // 14

group_cols = ['subject_id', 'two_week_period']
summary = data.groupby(group_cols).agg({'n_steps':'median', 
                                        'steps_ok':'mean'}) \
              .reset_index()

value_vars = ['n_steps', 'steps_ok']
summary = pd.melt(summary, id_vars=group_cols, value_vars=value_vars, 
                  var_name='feature_digital', value_name='value_digital')

summary['qc_type'] = 'steps'

summary.to_csv(output_file, index=False)
