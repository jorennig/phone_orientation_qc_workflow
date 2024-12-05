import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'qc_value': str}
data = pd.read_csv(input_file, dtype=dtype)

height = 5*(data['feature_digital'].nunique())
width = 5*(data['qc_type'].nunique())
plot = p9.ggplot(data, p9.aes(x='qc_value', y='value_digital')) \
     + p9.geom_boxplot(outlier_alpha=0.0) \
     + p9.geom_jitter(width=0.3, height=0.0, size=1.0, stroke=0.0, alpha=0.05) \
     + p9.theme(subplots_adjust={'wspace': 0.15}) \
     + p9.facet_grid('feature_digital ~ qc_type', scales='free_y') \
     + p9.xlab('') \
     + p9.ylab('')

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
