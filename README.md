# SpectralClustering
This repository provides functions that implement Spectral Clustering and Consensus Spectral Clustering in MatLab, and handle scripts for using the functions.

The function **spectralclustering.m** implements the spectral clustering algorithm derived from the classical algorithm published in 2002 by Ng. This function identifies clusters from data represented by an affinity matrix (rather than feature vectors as in the standard version). The affinity between each pair of data points is defined by the user.


The function takes in seven parameters, which include:

**affinityMatrix**: a n*n affinity matrix to represent the data.

**dataNameList**: a n*1 cell array to indicate the name (for e.g. gene name) of each data point.

**k**: the number of clusters to be identified.

**sigmag**: the standard deviation of Gaussian kernel to be used to transfer the input matrix.

**threshold**: the threshold to be used to filter the transformed affinities. Edges with transformed affinities below the thresholds would be removed.

**figureShow**: a logical value to indicate whether to show the figures for intermediate results.

**distanceTransform**: a integer to indicate how the input matrix should be transformed into a distance matrix (0 for no transformation, if the input matrix is already a distance matrix rather than an affinity matrix; 1 for taking inverse, 2 for taking the negative).


The function returns two variables:

**dataNameList**: a m*1 cell array providing the name of the remaining data after doing the filtering and clustering.

**clusterIndex**: a m*1 vector providing the cluster index of each remaining data point. The vector is arranged in the same order as the dataNameList.


The function **spectralclustering_consensus.m** implements the consensus extension of spectral clustering. It takes in one more parameter, **repeat**, an integer value indicates the number of replications of clustering. The consensus matrix will be visualized after calling the function.


**scUserHandle.m** and **scUserHandle_consensus.m** provide two example handles for using the functions and writing the results into a file under the Results directory. Users should change the input and output file directories to the ones on their local machine.

An example output of running scUserHandle.m would be

0 data points have been removed.

spectral clustering completed.

cluster1:431

cluster2:135

cluster3:223

cluster4:301

cluster5:214

cluster6:288

cluster7:155
