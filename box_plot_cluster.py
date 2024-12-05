import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
feature = snakemake.params.feature
value = snakemake.params.value
cluster_col = snakemake.params.cluster_col
alpha = snakemake.params.alpha
output_file = snakemake.output[0]

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)

plot = p9.ggplot(data, p9.aes(x=cluster_col, y=value)) \
     + p9.geom_boxplot(outlier_alpha=0.0) \
     + p9.geom_jitter(width=0.3, height=0.0, size=1.0, stroke=0.0, alpha=alpha) \
     + p9.theme(subplots_adjust={'wspace': 0.15}) \
     + p9.facet_wrap(feature, nrow=1, scales='free_y')

height = 5
width = height*(data[feature].nunique())

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
