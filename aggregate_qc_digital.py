import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

data['two_week_period'] = data['day_of_study'] // 14

cols = ['subject_id', 'two_week_period', 'feature_digital', 'qc_type']
aggregated = data.groupby(cols)['value_digital'].median() \
                 .reset_index()

aggregated.to_csv(output_file, index=False)
