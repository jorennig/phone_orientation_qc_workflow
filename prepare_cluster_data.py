import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

features = ['main_cluster', 'changes_percent_cluster']
data = data[data['feature_orientation'].isin(features)]

data = pd.pivot_table(data, index=['subject_id', 'two_week_period'], 
                      columns='feature_orientation', 
                      values='value_orientation') \
         .reset_index()

data.to_csv(output_file, index=False)
