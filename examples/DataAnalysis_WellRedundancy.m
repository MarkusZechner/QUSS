
% this script analyses the data (simulation output) to identify
% redundance in wells


close all; clear all; clc;
addpath('../src/DataAnalysis/');
addpath('../src/thirdparty/export_fig/');

prior_path = '../data/evidential/PriorData.mat';
load(prior_path);

% Saving Figures
BasePath = 'C:\Users\Markus Zechner\Documents\GitHub\QUSS\figures\DataAnalysis';
Case = 'OMV8TH_largerPrior';
% Trial is a subfolder of Case
Trial = 'Test1';

SaveFig = 'On';
[ FigureFolder ] = creatingFigureFolder(SaveFig, BasePath, Case, Trial);




%% Removing the observed data

TruthRealization = 1;
NumPriorRealizations=length(PriorData.data);
AvailableRealizations = setdiff(1:NumPriorRealizations,TruthRealization);

PriorData.data = PriorData.data(AvailableRealizations,:,:);

%% Performing PCA
% these values are only used to demonstrate the reconstruction of the
% signal - ALL components are ues for further analysis!!!
MinEigenValues = 3; EigenTolerance = 0.95;

rmpath('../src/thirdparty/fda_matlab/');
PlotLevel = 1; % 1=yes, 2=no
PriorDataPCA = ComputePCA(PriorData, EigenTolerance, MinEigenValues, PlotLevel, ...
    'Data', FigureFolder );


%% Performing CCA on between all wells

NbWells = size(PriorData.data, 3);
NbOperations = NbWells * NbWells;
k = 1;
for i=1:NbWells
    for j=1:NbWells
        
        
        d_f_well_i = PriorDataPCA{i}.score;
        d_f_well_j = PriorDataPCA{j}.score;
        
        fprintf('Working on well pair %d of %d .\n',k,NbOperations);
        
        [A, B, ~, d_c_well_i,d_c_well_j] = canoncorr(d_f_well_i,d_f_well_j);
        
        rho_temp = corrcoef(d_c_well_i(:,1),d_c_well_j(:,1));
        rho(i,j) = rho_temp(1,2);
        
        k = k + 1;
        
        
        
    end
end

%% PCA on the rho matrix
[coeff,score,latent,tsquared,explained,mu] = pca(rho);

cumexplained = cumsum(explained)/sum(explained);
plot(cumexplained)
xlabel('Number of Eigenvalues');
ylabel('Variance Explained');

scatter3(score(:,1),score(:,2),score(:,3))

%% Visulize the rho matrix
%range = [0 100]; 
imagesc(rho)
colorbar

%% getting the cdf for rho
UpperRho = triu(rho);

rho_vector = sort(reshape(UpperRho,1, []));

NonZeroIndex = find(rho_vector>0);
rho_vector = rho_vector(NonZeroIndex);
cdfplot(rho_vector)
xlabel('Correlation Coefficient of Canonical Correlation in 1st dimension')
size(rho_vector)

%testumnner = (NbOperations-1)/2

plot(cumsum(rho_vector));

hist(rho_vector)

HowManyAbove = size(find(rho_vector>0.9));



