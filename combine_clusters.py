import pandas as pd

input_file = snakemake.input.data
combine_clusters = snakemake.params.combine_clusters
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype = {'subject_id': str})

data['cluster'] = data['cluster'].replace(combine_clusters)

data.to_csv(output_file, index=None, header=True)
