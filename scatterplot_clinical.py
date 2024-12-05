import pandas as pd
from plotnine import ggplot, aes, geom_point, facet_wrap, xlab, ylab

input_file = snakemake.input.data[0]
input_results = snakemake.input.results[0]
feature = snakemake.wildcards.feature
feature_type = snakemake.params.feature_type
value_type = snakemake.params.value_type
merge_cols = snakemake.params.merge_cols
p_threshold = snakemake.params.p_threshold
output_file = snakemake.output[0]

dtype = {'subject_id': str}
data = pd.read_csv(input_file, dtype=dtype)
data = data[data[feature_type]==feature]

results = pd.read_csv(input_results)

results['rho_str'] = 'rho=' + results['rho'].round(2).astype(str)
results['p_str'] = 'p=' + results['p'].round(4).astype(str)
results['n_str'] = 'n=' + results['n'].astype(int).astype(str)
results.loc[results['p']<p_threshold, 'p_str'] = 'p<' + str(p_threshold)
results['header'] = results['qc_type'] + ', ' + results['feature_clinical'] + \
                    '\n ' + results['rho_str'] + ', ' + results['p_str'] + \
                    ', ' + results['n_str']

data = data.merge(results, on=merge_cols)

qc_types = data['qc_type'].nunique()
height = 5*qc_types
width = 5*data['feature_clinical'].nunique()
plot = ggplot(data, aes(x='value_clinical', y=value_type)) \
        + geom_point(alpha=0.20, stroke=0, size=3) \
        + facet_wrap('header', nrow=qc_types, scales='free_x') \
        + xlab('') \
        + ylab(feature) \

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
