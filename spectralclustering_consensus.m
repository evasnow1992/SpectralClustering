function [dataNameList,clusterIndex] = spectralclustering_consensus(affinityMatrix,dataNameList,k,sigmag,threshold,distanceTransform, repeat)
% This function implements consensus spectral clustering
% Parameter:
% affinityMatrix, the n*n affinity matrix for data points
% dataNameList, the n*1 cell array of data names
% k, the number of clusters
% sigmag, standard deviation of the gaussian kernel
% threshold, the threshold used to filter edges after kernel
% transformation
% distanceTransform, a integer value indicates how the input affinity
% matrix should be transformed into a distance matrix (0 for no
% transformation, 1 for taking inverse, 2 for taking minus)
% repeat, an integer value indicates the number of replications of
% clustering
% If any parameter is missing, default value will be used
% Return:
% dataNameList, a m*1 cell array providing the list of name of remaining
% data after clustering, m<=n
% clusterIndex, a m*1 vector providing the cluster index for each data point
if nargin < 7
    repeat = 10;
end
if nargin < 6
    distanceTransform = 1;
end
if nargin < 5
    threshold = 0.3;
end
if nargin < 4
    sigmag = 0.1;
end
if nargin < 3
    k = 5;
end
if nargin < 2
    filename = '../Results/matlabSpectralClustering/matlab_deg_list_GBM.csv';
    dataNameList = textread(filename,'%s');
end
if nargin == 0
    filename = '../Results/matlabSpectralClustering/matlab_deg_matrix_GBM.csv';
    affinityMatrix = csvread(filename,1,1);
end

nonZero = affinityMatrix > 0;

%**** Use Gussian kernel to transform the affinity matrix
%** Transform the affinity matrix to distance matrix
if distanceTransform == 1
    affinityMatrix(nonZero) = 1./affinityMatrix(nonZero); % Take the inverse of affinity as distance
elseif distanceTransform == 2
    affinityMatrix(nonZero) = max(affinityMatrix(nonZero))-affinityMatrix(nonZero);
end
%mean(affinityMatrix(nonZero))
%std(affinityMatrix(nonZero))


%** Kernel transformation
affinityMatrix(nonZero) = exp(-affinityMatrix(nonZero).^2./(2*sigmag^2));
%mean(affinityMatrix(nonZero))
%std(affinityMatrix(nonZero))
%mean(affinityMatrix(:))
%std(affinityMatrix(:))


%** Get rid of low affinities
lowAffinity = affinityMatrix < threshold;
%sum(lowAffinity(:)) % How many cells are below the threshold
affinityMatrix(lowAffinity) = 0;

%**** Get Laplacian
beforeDim = size(affinityMatrix,1);
sumByRow = sum(affinityMatrix,2);
zeroIndex = transpose(find(not(sumByRow)));
sumByRow(zeroIndex,:)=[];
affinityMatrix(zeroIndex,:)=[];
affinityMatrix(:,zeroIndex)=[];
dataNameList(zeroIndex,:)=[];
afterDim = size(affinityMatrix,1);
disp(strcat(num2str(beforeDim-afterDim), ' data points have been removed.'));
D = diag(sumByRow);
nonZeroD = D > 0;
D(nonZeroD) = D(nonZeroD).^(-1/2);
L = D*affinityMatrix*D;

%**** Get the k largest eigenvectors
%opts.tol = 1e-3;
%[X,E] = eigs(L,k,'lr',opts);
[X,E] = eigs(L,k);

%**** Normalize the eigenvector matrix by row
Y = normr(X);

%**** Consensus K mean
N = size(dataNameList, 1);
consensusMatrix = zeros(N, N);
for ii = 1:repeat
    disp(['Consensus clustering repeat: ', num2str(ii)])
    clusterIndex = kmeans(Y,k);
    for a = 1:N
        for b = 1:N
            if clusterIndex(a)==clusterIndex(b)
                consensusMatrix(a,b) = consensusMatrix(a,b) + 1;
            end
        end
    end
end

%****Visualize the consensus matrix
%***Reshape the consensus matrix according to cluster indices
accClusterSize = zeros(1, k);
newConsensusMatrix = zeros(N, N);
newMatrixIndex = ones(1, k);
newGeneMatrixIndex = zeros(1, N);
for a = 1:N
    accClusterSize(1, clusterIndex(a)) = accClusterSize(1, clusterIndex(a)) + 1;
end
for a = 2:k
    accClusterSize(1, a) = accClusterSize(1, a) + accClusterSize(1, a-1);
end
for a = 2:k
    newMatrixIndex(1, a) = newMatrixIndex(1, a) + accClusterSize(1, a-1);
end
for a = 1:N
    newGeneMatrixIndex(1, a) = newMatrixIndex(1, clusterIndex(a));
    newMatrixIndex(1, clusterIndex(a)) = newMatrixIndex(1, clusterIndex(a)) + 1;
    index1 = newGeneMatrixIndex(1, a);
    for b = 1:a
        index2 = newGeneMatrixIndex(1, b);
        newConsensusMatrix(index1, index2) = consensusMatrix(a,b);
        newConsensusMatrix(index2, index1) = consensusMatrix(b,a);
    end
end
I = mat2gray(newConsensusMatrix);
figure
imshow(I, [])
title('Consensus Matrix')

disp('spectral clustering completed.')
return