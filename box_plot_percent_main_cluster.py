import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

value = 'changes_percent_cluster'
cluster_col = 'main_cluster'
alpha = 0.2

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)

plot = p9.ggplot(data, p9.aes(x=cluster_col, y=value)) \
     + p9.geom_boxplot(outlier_alpha=0.0) \
     + p9.geom_jitter(width=0.3, height=0.0, size=1.0, stroke=0.0, alpha=alpha) \
     + p9.theme(subplots_adjust={'wspace': 0.15})

height = 5
width = 5

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
