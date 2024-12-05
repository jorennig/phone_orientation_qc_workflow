import pandas as pd
from sklearn.cluster import KMeans

input_file = snakemake.input.data
output_file = snakemake.output[0]
n_clust = int(snakemake.wildcards.n_clust)
kmeans_kwargs = snakemake.params.kmeans_args

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

features = data[['x', 'y', 'z']]

kmeans = KMeans(n_clusters=n_clust, **kmeans_kwargs)
kmeans.fit(features)
cluster_labels = kmeans.labels_

data['cluster'] = cluster_labels.tolist()
data['cluster'] = data['cluster'] + 1

data.to_csv(output_file, index=False)
