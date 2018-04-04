function [pca_scores] = ConcatinatePCAMatrix(FunctionalStruct,truth_real,...
    NumberOfEigenValues, Normalize)
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
    
    %     if Normalize = true
    %         % Normalize the PCA scores by the first singular value, which is the
    %         norm_score = FunctionalStruct{i}.score/sqrt(FunctionalStruct{i}.latent(1));
    %     end
    
    norm_score_temp = FunctionalStruct{i}.score;
    
    norm_score = norm_score_temp(:,1:NumberOfEigenValues);
    
    % Concanate the norm_score
    norm_scores = [norm_scores norm_score];
end




avail_real = setdiff(1:size(norm_scores,1),truth_real);

%pca_obs = mpca_scores(truth_real, :);
pca_scores = norm_scores(avail_real,:);


[rows, col] = size(pca_scores);
fprintf('There are %d number of wells.\n',num_wells);
fprintf('The concatinated matrix has dimensions %d x %d.\n',rows, col);

end






