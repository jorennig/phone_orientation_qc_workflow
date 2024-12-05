import pandas as pd
from sklearn.decomposition import PCA

input_file = snakemake.input.data[0]
output_file = snakemake.output.data
output_variance_explained = snakemake.output.variance_explained

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

vals = data[['x', 'y', 'z']].values
pca = PCA()
principal_components = pca.fit_transform(vals)

pc_cols = list(f'PC_{i}' for i in range(1, len(vals[0])+1))

variance_explained = pd.DataFrame(pca.explained_variance_ratio_*100, 
                                  index=pc_cols, 
                                  columns=['variance_explained']) \
                       .reset_index().rename(columns={'index': 'pc'})

principal_components = pd.DataFrame(principal_components, columns=pc_cols)
data = pd.concat([data, principal_components], axis=1)

data.to_csv(output_file, index=None, header=True)
variance_explained.to_csv(output_variance_explained, index=None, header=True)
