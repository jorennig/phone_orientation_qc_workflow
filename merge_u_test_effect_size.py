import pandas as pd

u_test_file = snakemake.input.u_test[0]
effect_size_file = snakemake.input.effect_size[0]
feature = snakemake.params.feature
output_file = snakemake.output[0]

u_test = pd.read_csv(u_test_file)
effect_size = pd.read_csv(effect_size_file)

merged = u_test.merge(effect_size, on=[feature, 'group1', 'group2'])

merged.to_csv(output_file, index=None, header=True)
