# -*- coding: utf-8 -*-
from __future__ import annotations
import os,sys,warnings
from typing import Optional

import numpy as np
import scipy.sparse as sps
from scipy.linalg import issymmetric
from sklearn.base import clone
from sklearn.neighbors import KNeighborsTransformer, kneighbors_graph
from scipy.spatial.distance import squareform

from datafold.pcfold.distance import BruteForceDist
from datafold.pcfold import PCManifoldKernel


def kneighbors_graph_from_transformer(X, knn_transformer=KNeighborsTransformer, mode="connectivity", include_self=True, **transformer_kwargs):
    """
    Calculate distance or connectivity with self generalized to any KNN transformer
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples
            and n_features is the number of features.

        knn_transformer : knn_transformer-like, [Default: KNeighborsTransformer]
            Either a fitted KNN transformer or an uninstantiated KNN transformer 
            with parameters specified **kwargs

        mode : str [Default: distance]
            Type of returned matrix: ‘connectivity’ will return the connectivity matrix with ones and zeros, 
            and ‘distance’ will return the distances between neighbors according to the given metric.

        include_self: bool [Default: True]
            Whether or not to mark each sample as the first nearest neighbor to itself. 
            If ‘auto’, then True is used for mode=’connectivity’ and False for mode=’distance’.
            

        transformer_kwargs: 
            Passed to knn_transformer if not instantiated
            
        Returns
        -------
        knn_graph : array-like, shape (n_samples, n_samples),

            scipy.sparse.csr_matrix

    """

    # mode checks
    assert mode in {"distance", "connectivity"}, "mode must be either 'distance' or 'connectivity'"

    # include_self checks
    if include_self == "auto":
        if mode == "distance":
            include_self = False
        else:
            include_self = True

    # If not instantiated, then instantiate it with **transformer_kwargs
    if not isinstance(knn_transformer, type):
        # Get params from model and add n_neighbors -= 1
        if include_self:
            assert not bool(transformer_kwargs), "Please provide uninstantiated `knn_transformer` or do not provide `transformer_kwargs`"
            warnings.warn("`include_self=True and n_neighbors=k` is equivalent to `include_self=False and n_neighbors=(k-1). Backend is creating a clone with n_neighbors=(k-1)")
            knn_transformer = clone(knn_transformer)
            n_neighbors = knn_transformer.get_params("n_neighbors")
            knn_transformer.set_params({"n_neighbors":n_neighbors - 1})
    else:
        try:
            n_neighbors = transformer_kwargs["n_neighbors"]
        except KeyError:
            raise Exception("Please provide `n_neighbors` as kwargs (https://docs.python.org/3/glossary.html#term-argument)")
        if include_self:
            transformer_kwargs["n_neighbors"] = n_neighbors - 1
        knn_transformer = knn_transformer(**transformer_kwargs)
        
    # Compute KNN graph for distances
    knn_graph = knn_transformer.fit_transform(X)
    
    # Convert to connectivity
    if mode == "connectivity":
        # Get all connectivities
        knn_graph = (knn_graph > 0).astype(float)
           
        # Set diagonal to 1.0
        if include_self:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                knn_graph.setdiag(1.0)
                
    return knn_graph

def brute_force_kneighbors_graph_from_rectangular_distance(distance_matrix, n_neighbors:int, mode="connectivity", include_self=True):
    assert mode in {"distance", "connectivity"}, "mode must be either 'distance' or 'connectivity'"

    # ==================================================================
    # Naive
    # -----
    # out = csr_matrix(distance_matrix.shape)
    # if mode == "connectivity":
    #     for i, index in enumerate(indices):
    #         out[i,index] = 1.0
    # if mode == "distance":
    #     distances = np.sort(distance_matrix, axis=1)[:, :n_neighbors]
    #     for i, index in enumerate(indices):
    #         out[i,index] = distances[i]
    # ==================================================================
    if include_self:
        n_neighbors = n_neighbors - 1
        
    # Sort indices up to n_neighbors            
    indices = np.argpartition(distance_matrix, n_neighbors, axis=1)[:, :n_neighbors]
    # Use ones for connectivity values
    if mode == "connectivity":
        data = np.ones(distance_matrix.shape[0] * n_neighbors, dtype=float)
    # Use distances values
    if mode == "distance":
        data = np.partition(distance_matrix, n_neighbors, axis=1)[:, :n_neighbors].ravel()
    # Get row indices
    row = np.repeat(np.arange(distance_matrix.shape[0]), n_neighbors)
    # Get column indicies
    col = indices.ravel()
    
    # Build COO matrix
    graph = sps.coo_matrix((data, (row, col)), shape=distance_matrix.shape)

    # Convert to CRS matrix
    return graph.tocsr()
    
    
class KNeighborsKernel(PCManifoldKernel):
    """
    K-Nearest Neighbors Kernel
    
    Acknowledgement: 
    https://gitlab.com/datafold-dev/datafold/-/issues/166
    """
    def __init__(self, metric:str, n_neighbors:int, distance_matrix:Optional[np.ndarray]=None, copy_distance_matrix=False, verbose=0):

        self.n_neighbors = n_neighbors
        self.verbose = verbose
        self.copy_distance_matrix = copy_distance_matrix
        if distance_matrix is not None:
            if len(distance_matrix.shape) == 1:
                distance_matrix = squareform(distance_matrix)
            else:
                if copy_distance_matrix:
                    distance_matrix = distance_matrix.copy()
        self.distance_matrix = distance_matrix

        distance = BruteForceDist(metric=metric)
        super().__init__(is_symmetric=True, is_stochastic=False, distance=distance)

    def __call__(self, X: np.ndarray, Y: Optional[np.ndarray] = None, **kernel_kwargs):
        if all([
            Y is None,
            self.distance_matrix is not None,
            ]):
            n, m = X.shape
            assert self.distance_matrix.shape[0] == n, "X.shape[0] must equal distance_matrx.shape[0]"
            assert self.distance_matrix.shape[1] == n, "X.shape[0] must equal distance_matrx.shape[1]"
            distance_matrix = self.distance_matrix
            if self.verbose > 0:
                print("Precomputed distance matrix detected. Skipping pairwise distance calculations.", file=sys.stderr, flush=True)
        else:
            distance_matrix = self.distance(X, Y)
        return self.evaluate(distance_matrix)

    def evaluate(self, distance_matrix):

        # Compute KNN connectivity kernel
        distance_matrix_is_square = False
        shape = distance_matrix.shape
        if shape[0] == shape[1]:
            if issymmetric(distance_matrix):
                distance_matrix_is_square = True
        if distance_matrix_is_square:
            connectivities = kneighbors_graph(distance_matrix, n_neighbors=self.n_neighbors, metric="precomputed", include_self=True, mode="connectivity")
        else:
            connectivities = brute_force_kneighbors_graph_from_rectangular_distance(distance_matrix, n_neighbors=self.n_neighbors, include_self=True, mode="connectivity")

        return connectivities