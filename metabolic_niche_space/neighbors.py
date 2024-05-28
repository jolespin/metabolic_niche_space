# -*- coding: utf-8 -*-
from __future__ import annotations
import os,sys,warnings
from typing import Optional

import numpy as np
import scipy.sparse as sps
from scipy.linalg import issymmetric
from sklearn.base import clone
from sklearn.neighbors import KNeighborsTransformer, kneighbors_graph
# from scipy.spatial.distance import pdist, squareform

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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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
    def __init__(self, metric:str, n_neighbors:int):

        self.n_neighbors = n_neighbors
        distance = BruteForceDist(metric=metric)
        super().__init__(is_symmetric=True, is_stochastic=False, distance=distance)

    def __call__(self, X: np.ndarray, Y: Optional[np.ndarray] = None, **kernel_kwargs):
        distance = self.distance(X, Y)
        return self.evaluate(distance)

    def evaluate(self, distance_matrix):

        def _brute_force_knn_connectivity_from_rectangular_distance(distance_matrix, n_neighbors:int, include_self=True,  mode="connectivity"):
            assert mode in {"distance", "connectivity"}, "mode must be either 'distance' or 'connectivity'"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if include_self:
                    n_neighbors = n_neighbors - 1
                indices = np.argsort(distance_matrix, axis=1)[:, :n_neighbors]
        
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
        
                # Flatten the inds and dists arrays
                col_indices = indices.flatten()
        
                if mode == "connectivity":
                    data = np.ones_like(col_indices, dtype=float)
                if mode == "distances":
                    distances = np.sort(distance_matrix, axis=1)[:, :n_neighbors]
                    data = distances.flatten()
                
                # Create the row indices
                number_of_rows = indices.shape[0]
                row_indices = np.repeat(np.arange(number_of_rows), indices.shape[1])
                
                # Create the CSR matrix
                out = sps.csr_matrix((data, (row_indices, col_indices)), shape=distance_matrix.shape)
                
                return out

        # Compute KNN connectivity kernel
        distance_matrix_is_rectangular = True
        shape = distance_matrix.shape
        if shape[0] == shape[1]:
            if issymmetric(distance_matrix):
                distance_matrix_is_rectangular = False
        if distance_matrix_is_rectangular:
            connectivities = brute_force_kneighbors_graph_from_rectangular_distance(distance_matrix, n_neighbors=self.n_neighbors, include_self=True, mode="connectivity")
        else:
            connectivities = kneighbors_graph(distance_matrix, n_neighbors=self.n_neighbors, metric="precomputed", include_self=True, mode="connectivity")

            
        return connectivities