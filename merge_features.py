import pandas as pd

input_data = snakemake.input.data[0]
merge_cols = snakemake.params.merge_cols
input_test = snakemake.input.test[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_data, dtype={'subject_id': str})
test = pd.read_csv(input_test, dtype={'subject_id': str})

digital_cols = ['feature_digital', 'value_digital']
test = test.filter(merge_cols+digital_cols)

merged = data.merge(test, on=merge_cols)

merged.to_csv(output_file, index=None, header=True)
