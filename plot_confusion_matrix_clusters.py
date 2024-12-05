import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

data = pd.read_csv(input_file)

data = pd.melt(data, id_vars='main_orientation', value_vars=['1', '2', '3'],
               var_name='cluster', value_name='n')

order_clusters = ['1', '2', '3']
data['cluster'] = pd.Categorical(data['cluster'], ordered=True, 
                                 categories=order_clusters)

p = p9.ggplot(data, p9.aes(x='cluster', y='main_orientation', fill='n')) \
    + p9.geom_tile() \
    + p9.geom_text(p9.aes(label='n'), alpha=1, size=10) \
    + p9.scale_fill_gradient(low='white', high='red') \
    + p9.theme(legend_position='none')

p.save(output_file, dpi=300, width=10, height=10, verbose=False)
