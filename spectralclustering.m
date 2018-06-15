function [dataNameList,clusterIndex] = spectralclustering(affinityMatrix,dataNameList,k,sigmag,threshold,figureShow,distanceTransform)
% This function implements spectral clustering
% Parameter:
% affinityMatrix, the n*n affinity matrix for data points
% dataNameList, the n*1 cell array of data names
% k, the number of clusters
% sigmag, standard deviation of the gaussian kernel
% threshold, the threshold used to filter edges after kernel
% transformation
% figureShow, a logical value indicates whether to show the figures for
% intermediate results
% distanceTransform, a integer value indicates how the input affinity
% matrix should be transformed into a distance matrix (0 for no
% transformation, 1 for taking inverse, 2 for taking minus)
% If any parameter is missing, default value will be used
% Return:
% dataNameList, a m*1 cell array providing the list of name of remaining
% data after clustering, m<=n
% clusterIndex, a m*1 vector providing the cluster index for each data point
if nargin < 7
    distanceTransform = 1;
end
if nargin < 6
    figureShow = true;
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
%**** Visualize the initial affinity matrix (optional)
if figureShow == 1
    I = mat2gray(affinityMatrix);
    imshow(I,[])
    title('Heatmap for initial affinity matrix')
    %print('../Results/initialHeatmap','-djpeg','-r600')
    
    figure
    histogram(affinityMatrix(:),40)
    title('Histogram for initial affinity matrix')
    figure
    histogram(affinityMatrix(nonZero),40)
    title('Histogram for initial affinity matrix nonzero cells')
end

%**** Use Gussian kernel to transform the affinity matrix
%** Transform the affinity matrix to distance matrix
if distanceTransform == 1
    affinityMatrix(nonZero) = 1./affinityMatrix(nonZero); % Take the inverse of affinity as distance
elseif distanceTransform == 2
    affinityMatrix(nonZero) = max(affinityMatrix(nonZero))-affinityMatrix(nonZero);
end
%mean(affinityMatrix(nonZero))
%std(affinityMatrix(nonZero))

%** Visualize the distance matrix nonzero cells (optional)
if figureShow == 1
    figure
    yyaxis left
    histogram(affinityMatrix(nonZero),40)
    title('Histogram for distance matrix nonzero cells')
    hold on
    x1 = linspace(min(affinityMatrix(nonZero)), max(affinityMatrix(nonZero)));
    y1 = exp(-x1.^2./(2*sigmag^2));
    yyaxis right
    plot(x1, y1)
    hold off
    %print('../Results/distanceNonzero','-djpeg','-r600')
    
end

%** Kernel transformation
affinityMatrix(nonZero) = exp(-affinityMatrix(nonZero).^2./(2*sigmag^2));
%mean(affinityMatrix(nonZero))
%std(affinityMatrix(nonZero))
%mean(affinityMatrix(:))
%std(affinityMatrix(:))

%** Visualize the transformed affinity matrix (optional)
if figureShow == 1
    figure
    histogram(affinityMatrix(:),40)
    title('Histogram for affinity matrix after kernal transformation')
    hold on
    %plot([threshold, threshold],[0, 18*10^5])
    %print('../Results/affinityKernelTransfer','-djpeg','-r600')
    
    figure
    histogram(affinityMatrix(nonZero),40)
    title('Histogram for affinity matrix nonzero cells after kernal transformation')
end

%** Get rid of low affinities
lowAffinity = affinityMatrix < threshold;
%sum(lowAffinity(:)) % How many cells are below the threshold
affinityMatrix(lowAffinity) = 0;

%** Visualize the filtered (final) affinity matrix (optional)
if figureShow == 1
    I = mat2gray(affinityMatrix);
    figure
    imshow(I, [])
    title('Heatmap for final affinity matrix')
    %print('../Results/finalHeatmap','-djpeg','-r600')
    
    figure
    histogram(affinityMatrix(:),40)
    title('Histogram for final affinity matrix')
    nonZero2 = affinityMatrix > 0;
    figure
    histogram(affinityMatrix(nonZero2),40)
    title('Histogram for final affinity matrix nonzero cells')
end

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

%**** K mean
clusterIndex = kmeans(Y,k);

disp('spectral clustering completed.')
return