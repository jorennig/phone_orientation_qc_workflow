import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
feature_test = snakemake.wildcards.feature_test
feature = snakemake.params.feature
value = snakemake.params.value
cluster_col = snakemake.params.cluster_col
output_file = snakemake.output[0]

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)
data = data[data[feature]==feature_test]

height = 5
width = height*(data[cluster_col].nunique())
plot = p9.ggplot(data, p9.aes(x=value)) \
         + p9.geom_histogram() \
         + p9.facet_wrap(cluster_col) \
         + p9.xlab(feature_test) \
         + p9.ggtitle(cluster_col)

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
