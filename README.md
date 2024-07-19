# Metabolic Niche Space

Coming soon...

#### Test KNN Kernel:

```python
# Built-In
import warnings

# PyData
from scipy.spatial.distance import pdist, squareform
from sklearn.datasets import make_classification
from sklearn.neighbors import kneighbors_graph, KNeighborsTransformer
from sklearn.model_selection import train_test_split

# Metabolic Niche Space
from metabolic_niche_space.neighbors import KNeighborsKernel
from metabolic_niche_space.manifold import DiffusionMaps # Shortcut: from datafold.dynfold import DiffusionMaps

# Create dataset
n_samples = 1000
n_neighbors=int(np.sqrt(n_samples))
X, y = make_classification(n_samples=n_samples, n_features=10, n_classes=2, n_clusters_per_class=1, random_state=0)
X_training, X_testing, y_training, y_testing = train_test_split(X, y, test_size=0.3)

# Build KNeighbors Kernel
kernel = KNeighborsKernel(metric="euclidean", n_neighbors=30)

# Calculate Diffusion Maps using KNeighbors
model = DiffusionMaps(kernel=kernel, n_eigenpairs=int(np.sqrt(X_training.shape[0])))
dmap_X = model.fit_transform(X_training)
dmap_Y = model.transform(X_testing)

# Shapes
print(dmap_X.shape, dmap_Y.shape)
# (700, 26) (300, 26)

# Plot
fig, ax = plt.subplots(figsize=(5,5))
ax.scatter(dmap_X[:,1], dmap_X[:,2], c=list(map(lambda c: {True:"red", False:"black"}[c], y_training)), alpha=0.1, label="Training")
ax.scatter(dmap_Y[:,1], dmap_Y[:,2], c=list(map(lambda c: {True:"red", False:"black"}[c], y_testing)), s = 5, alpha=0.75, label="Testing")
ax.legend()
```

![](images/test_knn.png)