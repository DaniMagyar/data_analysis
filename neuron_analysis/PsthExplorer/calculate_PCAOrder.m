function [leafOrder, Dend, Clusters] = calculate_PCAOrder(PSTHall, preferences)

% PCA analysis ------------------------------------------------------------
% In "PSTHall" the rows are the neurons
[~,PCA1]  = pca(PSTHall); % running pricipal component analysis
PCA2      = PCA1(:,1:preferences.PCA_num); % extracting the first PCA_num pricinpal components 
inNAN = nanmedian(PCA2,1); % to avoid NaN values, set them to the median
for jj = 1: size(PSTHall,1)
    if isnan(PCA2(jj,1)) == 1
        PCA2(jj,1:3) = inNAN;
    end
end
Dend      = linkage(PCA2,preferences.LKmethod,preferences.LKmetric); % calculating the dendrogram
Clusters  = cluster(Dend,'maxclust',preferences.clustnum); % clustering based on the tree
D         = pdist(PCA2); %euclidean distrance between point in the artifical space
leafOrder = optimalleaforder(Dend,D); % optimal order for plotting 
PCA        = PCA2;
Dendrogram = Dend;
