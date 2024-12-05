import pandas as pd
import plotnine as p9

input_file = snakemake.input.data[0]
qc_check_new = snakemake.wildcards.qc_check_new
qc_check_old = snakemake.wildcards.qc_check_old
output_file = snakemake.output[0]

dtype = {qc_check_old: str, qc_check_new: str}
data = pd.read_csv(input_file, dtype=dtype)

p = p9.ggplot(data, p9.aes(x=qc_check_new, y=qc_check_old, fill='n')) \
    + p9.geom_tile() \
    + p9.geom_text(p9.aes(label='n'), alpha=1, size=10) \
    + p9.scale_fill_gradient(low='white', high='red') \
    + p9.theme(legend_position='none')

p.save(output_file, dpi=300, width=10, height=10, verbose=False)
