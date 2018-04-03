function [mpca_scores, mpca_obs] = MixedPCAnoBsplines(FunctionalStruct,truth_real,...
    eigentolerance, SavePath)
%MIXEDPCA Computes Mixed PCA of a PCA Object
%
% Inputs:
%   FunctionalStruct: Structure containing harmscr for each response
%   variable
%   truth_real[Optional]: Whether to set aside a realization as d_obs
%   eigentolerance: Number of eigenvalues to keep to maintain this variance
% Outputs:
%   mpca_scores: Scores of response variables with 99% of variance kept
%   mpca_obs: Score of observed data
%
% Author: Lewis Li (lewisli@stanford.edu)
% updated for no B-splines Markus Zechner (mzechner@stanford.edu)
% Date: May 24th 2016
% Update: March 24th 2018

if (nargin < 4)
    SaveOn = false;
else
    SaveOn = true;
end

% This is the number of variables
num_wells = length(FunctionalStruct);

% Concenated normalized scores
norm_scores = [];

for i = 1:num_wells
    % Perform regular PCA on each well
    %[~,~,latent] = pca(FunctionalStruct{i}.score);
    
    % Normalize the PCA scores by the first singular value, which is the
    norm_score = FunctionalStruct{i}.score/sqrt(FunctionalStruct{i}.latent(1));
    
    % Concanate the norm_score
    norm_scores = [norm_scores norm_score];
end

% Perform PCA on concatenated matrix
[~,mpca_scores,~,~,explained] = pca(norm_scores);

% Compute explained variance
explained = cumsum(explained)/sum(explained);


plot(explained,'LineWidth',3)
set(gcf,'color','w');
%ylim([0 1])
xlabel('Number of Eigencomponents');
ylabel('Variance Explained');

if SaveOn == true
    export_fig([SavePath 'VarianceExplainedMixedPCA'], '-png','-m3');
end

% Check number of components to keep
eigenToKeep = 3;

if (eigentolerance<1)
    ix = max(find(explained > eigentolerance, 1, 'first'),eigenToKeep);
else
    ix = eigentolerance;
end

fprintf('%d Eigenvalues are use to describe %.2f of variance.\n',ix,eigentolerance);

% Whether we set aside a truth realization
if truth_real==0
    mpca_scores = mpca_scores(:,1:ix);
    mpca_obs = 0;
else
    avail_real = setdiff(1:size(mpca_scores,1),truth_real);
    mpca_obs = mpca_scores(truth_real,1:ix);
    mpca_scores = mpca_scores(avail_real,1:ix);
end

end

