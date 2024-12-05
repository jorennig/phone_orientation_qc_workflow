import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
dimension = snakemake.wildcards.dimension
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

height = 5
width = 5
plot = p9.ggplot(data, p9.aes(x=dimension)) \
         + p9.geom_histogram(bins=100)

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
