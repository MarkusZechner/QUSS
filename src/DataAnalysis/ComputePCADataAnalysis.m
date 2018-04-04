function [ PriorDataPCA ] = ComputePCADataAnalysis(PriorData, eigentolerance,eigenToKeep, PlotLevel, ...
    Type, FigureFolder)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Author: Markus Zechner (mzechner@stanford.edu)
% Date: March 24th 2018

if (nargin < 6)
    SaveOn = false;
else
    SaveOn = true;
end

PriorDataPCA = {};

% Number of wells
num_wells = length(PriorData.ObjNames);

NbSim = size(PriorData.data,1);


for i = 1:num_wells
    % Perform regular PCA on each well
    
    [coeff,score,latent,tsquared,explained,mu] = pca(PriorData.data(:,:,i));
    
    PriorDataPCA{i}.coeff = coeff;
    PriorDataPCA{i}.score = score;
    PriorDataPCA{i}.latent = latent;
    PriorDataPCA{i}.explained = explained;
    PriorDataPCA{i}.mu = mu;
    
end


if PlotLevel == 1
    
    RandomWell = randi([1 num_wells],1);
    RandomRealization = randi([1 NbSim],1);
    
    % Plotting the variance explained
    figure;
    explainedWell = cumsum(PriorDataPCA{RandomWell}.explained)/sum(PriorDataPCA{RandomWell}.explained);
    plot(explainedWell,'LineWidth',3)
    set(gcf,'color','w');
    %ylim([0 1])
    xlabel('Number of Eigencomponents');
    ylabel('Variance Explained');
    
    if SaveOn == true
        figure_name=['VarianceExplained_' Type];
        export_fig([FigureFolder figure_name], '-png','-m3');
    end
    
    % Computing how many Eigenvalues are needed for a given variance
    % explained
    
    if (eigentolerance<1)
        ix = max(find(explainedWell > eigentolerance, 1, 'first'),eigenToKeep);
    else
        ix = eigentolerance;
    end
    
    % undoing PCA
    
    Xhat = PriorDataPCA{RandomWell}.score(:,1:ix) * PriorDataPCA{RandomWell}.coeff(:,1:ix)';
    Xhat = bsxfun(@plus, Xhat, PriorDataPCA{RandomWell}.mu);
    
    
    figure;
    plot(PriorData.time, PriorData.data(RandomRealization,:,RandomWell), ...
        'o','MarkerFaceColor','b','MarkerEdgeColor','b');
    hold on;
    plot(PriorData.time, Xhat(RandomRealization,:), ...
        'o','MarkerFaceColor','r','MarkerEdgeColor','r');
    title({[ 'Reconstruction of Well ' PriorData.ObjNames{RandomWell} 'in Realization ', ...
        num2str(RandomRealization)]; [' with ' num2str(eigentolerance) ' explained variance']} , ...
        'FontSize',10)
    set(gcf,'color','w');
    xlabel('Time (days)');
    ylabel(PriorData.name);
    
    if SaveOn == true
        figure_name=['Reconstruction' Type];
        export_fig([FigureFolder figure_name], '-png','-m3');
    end
    
end

end

