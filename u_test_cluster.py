import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.stats.multicomp as mc

def u_tests_pairwise(chunk):    
    comp = mc.MultiComparison(chunk[value], chunk[cluster_col])
    tbl, a1, a2 = comp.allpairtest(mannwhitneyu, method='bonf')
    results_as_html = tbl.as_html()
    results = pd.read_html(results_as_html, header=0)[0]
    return results.filter(['group1', 'group2', 'stat', 'pval'])

input_file = snakemake.input.data[0]
feature = snakemake.params.feature
value = snakemake.params.value
cluster_col = snakemake.params.cluster_col
output_file = snakemake.output[0]

dtype = {'subject_id': str, cluster_col: str}
data = pd.read_csv(input_file, dtype=dtype)

u_test = data.groupby([feature]).apply(u_tests_pairwise) \
             .reset_index().drop(columns=['level_1'])

u_test.to_csv(output_file, index = None, header = True)
