import pandas as pd

input_file = snakemake.input.data[0]
qc_check_new = snakemake.wildcards.qc_check_new
qc_check_old = snakemake.wildcards.qc_check_old
output_file = snakemake.output[0]

dtype = {qc_check_old: str, qc_check_new: str}
data = pd.read_csv(input_file, dtype=dtype)

table = data[[qc_check_old, qc_check_new]].value_counts() \
                                          .reset_index(name='n')

table.to_csv(output_file, index=True, header=True)
