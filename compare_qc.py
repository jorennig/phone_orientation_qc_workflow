import pandas as pd
import numpy as np

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str}
data = pd.read_csv(input_file, dtype=dtype)

index_cols = ['subject_id', 'two_week_period', 'feature_digital']
data = pd.pivot_table(data, index=index_cols, columns='qc_type', 
                      values='value_digital') \
         .reset_index() \
         .dropna()

data['diff_value'] = data['phone_not_on_table_ok'] - data['new_qc_pass']
data['diff_value_log'] = np.log10(data['diff_value'])
data['diff_digital'] = data['diff_value'] != 0

data.to_csv(output_file, index=False)
