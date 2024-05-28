
# from .neighbors import KNeighborsKernel
from datafold.dynfold import DiffusionMaps, Roseland

class DiffusionMapEmbedding(DiffusionMaps):
    """
    DiffusionMapEmbedding is a renamed version of the DiffusionMaps class from the datafold.dynfold package.
    It inherits all methods and properties of the original class without any modifications.
    
    Documentation: 
        https://datafold-dev.gitlab.io/datafold/api/datafold.dynfold.DiffusionMaps.html
    Citation: 
        Lehmberg et al., (2020). datafold: data-driven models for point clouds and time series on manifolds. 
        Journal of Open Source Software, 5(51), 2283, https://doi.org/10.21105/joss.02283
    """
    pass

class LandmarkDiffusionMapEmbedding(Roseland):
    """
    LandmarkDiffusionMapEmbedding is a renamed version of the Roseland class from the datafold.dynfold package.
    It inherits all methods and properties of the original class without any modifications.
    
    Documentation: 
        https://datafold-dev.gitlab.io/datafold/api/datafold.dynfold.Roseland.html
    Citation: 
        Lehmberg et al., (2020). datafold: data-driven models for point clouds and time series on manifolds. 
        Journal of Open Source Software, 5(51), 2283, https://doi.org/10.21105/joss.02283
    """
    pass


__all__ = ["DiffusionMapEmbedding", "LandmarkDiffusionMapEmbedding"]