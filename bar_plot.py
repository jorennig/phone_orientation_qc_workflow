import pandas as pd
import plotnine as p9

input_file = snakemake.input.data
x = snakemake.params.x
y = snakemake.params.y
output_file = snakemake.output[0]

data = pd.read_csv(input_file)

height = 5
width = 5
plot = p9.ggplot(data, p9.aes(x=x, y=y)) \
         + p9.geom_bar(stat='identity') \

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
