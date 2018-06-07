clear
close all

filename = '../Results/Results_GBM/gene_list_GBM.csv';
degList = textread(filename,'%s');
filename = '../Results/Results_GBM/gene_matrix_GBM.csv';
degMatrix = csvread(filename,1,1);
k = 15;
sigmag = 0.1;
threshold = 0;
figureShow = true;
distanceTransform = 1;

clusterSize = zeros(1, k);
[degList, clusterIndex] = spectralclustering(degMatrix,degList,k,sigmag,threshold,figureShow, distanceTransform);

outfilename = strcat('../Results/spectralClustering_',num2str(k),'clusters_GBM.csv');
fid = fopen(outfilename, 'w');
for a = 1:size(degList,1)
    fprintf(fid, '%s,%d\n', degList{a},clusterIndex(a));
    clusterSize(1, clusterIndex(a)) = clusterSize(1, clusterIndex(a)) + 1;
end
fclose(fid) ;
for a = 1:k
    disp(strcat('cluster', int2str(a), ': ', int2str(clusterSize(1, a))));
end