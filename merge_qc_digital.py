import pandas as pd

input_qc = snakemake.input.qc[0]
input_digital = snakemake.input.digital[0]
merge_cols = snakemake.params.merge_cols
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str}
qc = pd.read_csv(input_qc, dtype=dtype)

cols = merge_cols + ['feature_digital', 'value_digital']
digital = pd.read_csv(input_digital, dtype=dtype) \
            .filter(cols)

merged = qc.merge(digital, on=merge_cols)

merged = pd.melt(merged, id_vars=cols, value_vars=
                 ['phone_not_on_table_ok', 'steps_ok', 'new_qc_pass'],
                 var_name='qc_type', value_name='qc_value')

merged.to_csv(output_file, index=None, header=True)
