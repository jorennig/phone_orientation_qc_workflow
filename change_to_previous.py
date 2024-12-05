import pandas as pd

input_file = snakemake.input.data[0]
variable = snakemake.params.variable
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

cols = ['subject_id']
data['shift'] = data.groupby(cols)[variable].shift(1)
data['shift'] = data['shift'].fillna(method='bfill')

data['change_to_previous'] = (data[variable] != data['shift']).astype(int)

data['changes_total_' + variable] = data.groupby(cols)['change_to_previous'].transform('sum')

data['changes_percent_' + variable] = data['changes_total_' + variable] / \
                                      data['n_tests'] * 100

data = data.drop(columns=['shift'])

data.to_csv(output_file, index=None, header=True)
