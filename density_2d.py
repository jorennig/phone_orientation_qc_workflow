import pandas as pd
from plotnine import ggplot, aes, geom_point, geom_density_2d, xlim, ylim

input_file = snakemake.input.data[0]
coords = snakemake.wildcards.coords
c = snakemake.params.coordinates_dict[coords]
color = snakemake.params.color
output_file = snakemake.output[0]

dtype = {'subject_id': str, 'device_id': str, 'cluster': str}
data = pd.read_csv(input_file, dtype=dtype)

width = 5
height = 5
plot = ggplot(data, aes(x=c[0], y=c[1], color=color)) \
        + geom_point(alpha=0.002, stroke=0, size=2) \
        + geom_density_2d() \
        + xlim(0, 10) \
        + ylim(0, 10)

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
