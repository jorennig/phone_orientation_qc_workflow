import pandas as pd
from sklearn.cluster import KMeans
from kneed import KneeLocator
from plotnine import ggplot, aes, geom_line, annotate

input_file = snakemake.input.data
output_file = snakemake.output[0]
kmeans_kwargs = snakemake.params.kmeans_args
max_clust = int(snakemake.params.max_clust)

dtype = {'subject_id': str, 'device_id': str}
data = pd.read_csv(input_file, dtype=dtype)

features = data[['x', 'y', 'z']]

sse = []
for k in range(1, max_clust):
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(features)
    sse.append(kmeans.inertia_)

n_clusters = (list(range(1, max_clust)))

data = pd.DataFrame({'n_clusters': n_clusters, 'SSE': sse})
data['n_clusters'] = pd.Categorical(data['n_clusters'], ordered=True)

kl = KneeLocator(range(1, max_clust), sse, curve='convex', direction='decreasing')
e = kl.elbow
s = kl.elbow_y

width = 5
height = 5
plot = ggplot(data, aes(x='n_clusters', y='SSE', group=1)) \
        + geom_line() \
        + annotate('point', x=e, y=s, colour = 'red')

plot.save(output_file, dpi=300, width=width, height=height, verbose=False)
