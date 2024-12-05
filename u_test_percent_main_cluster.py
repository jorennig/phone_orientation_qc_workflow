import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.stats.multicomp as mc

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]

value = 'changes_percent_cluster'
cluster_col = 'main_cluster'

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)

comp = mc.MultiComparison(data[value], data[cluster_col])
tbl, a1, a2 = comp.allpairtest(mannwhitneyu, method='bonf')
results_as_html = tbl.as_html()
u_test = pd.read_html(results_as_html, header=0)[0].filter(['group1', 'group2', 'stat', 'pval'])

u_test.to_csv(output_file, index = None, header = True)
