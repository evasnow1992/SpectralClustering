# SpectralClustering
The function spectralclustering.m implements the spectral clustering algorithm derived from the classical algorithm published in 2002 by Ng.


The function takes in seven parameters, which include:

affinityMatrix: a n*n matrix to represent the graph.

dataNameList: a n*1 cell array to indicate the name (for e.g. gene name) of the data.

k: the number of clusters.

sigmag: the standard deviation to be used in the gaussian kernel in transforming the input matrix.

threshold: the threshold to be used in filtering transformed affinities. Edges with transformed affinities below the thresholds would be removed.

figureShow: a logical value to indicate whether to show the figures for intermediate results.

distanceTransform: a integer to indicate how the input matrix should be transformed into a distance matrix (0 for no transformation, used when the input matrix is already a distance matrix rather than an affinity matrix; 1 for taking inverse, 2 for taking the negative).


The function returns two variables:

dataNameList: a m*1 cell array providing the name of the remaining data after doing the filtering and clustering.

clusterIndex: a m*1 vector providing the cluster index of each remaining data point. The vector is arranged in the same order as the dataNameList.



scUserHandle_BRCA_metabric.m and scUserHandle_GBM.m provide two examples in using the function spectralclustering and writing the results into a file under the Results directory. Users should change the input and output file directories to the ones on their local machine.

An example output of running scUserHandle_BRCA_metabric.m would be

0 data points have been removed.

spectral clustering completed.

cluster1:431

cluster2:135

cluster3:223

cluster4:301

cluster5:214

cluster6:288

cluster7:155
