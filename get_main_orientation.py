import pandas as pd

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

cols = ['subject_id']
data['n_tests'] = data.groupby(cols)['subject_id'].transform('count')

data['main_orientation'] = data[['x', 'y', 'z']].apply(lambda x: x.argmax(), axis=1)
data['main_orientation'] = data['main_orientation'].replace(dict(enumerate('xyz')))

data['main_orientation_total'] = data.groupby(cols)['main_orientation'].transform(lambda x: x.value_counts().idxmax())
data['main_orientation_total_n'] = data.groupby(cols)['main_orientation'].transform(lambda x: x.value_counts().max())

data['percent_main_orientation_total'] = data['main_orientation_total_n'] / \
                                         data['n_tests'] * 100

data.to_csv(output_file, index=None, header=True)
