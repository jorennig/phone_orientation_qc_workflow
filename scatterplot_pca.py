import pandas as pd
from plotnine import ggplot, aes, geom_point

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str, 'cluster': str}
data = pd.read_csv(input_file, dtype=dtype)

width = 5
height = 5
plot = ggplot(data, aes(x='PC_1', y='PC_2', color='cluster')) \
        + geom_point(alpha=0.002, stroke=0, size=2)

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
