import pandas as pd

input_file = snakemake.input.data[0]
features = snakemake.params.features
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

data = data.filter(['subject_id', 'n_tests'] + features) \
           .drop_duplicates()

data = pd.melt(data, id_vars=['subject_id', 'n_tests'], 
               value_vars=features,
               var_name='feature_orientation', value_name='value_orientation')

data['two_week_period'] = 0
data['qc_type'] = 'none'

data.to_csv(output_file, index=None, header=True)
