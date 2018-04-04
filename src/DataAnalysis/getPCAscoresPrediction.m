function [pca_scores_well] = getPCAscoresPrediction(FunctionalStruct,truth_real,...
    NumberOfEigenValues,WellNumber)
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


norm_score_temp = FunctionalStruct{WellNumber}.score;

PredictionScoreWell = norm_score_temp(:,1:NumberOfEigenValues);


avail_real = setdiff(1:size(PredictionScoreWell,1),truth_real);

%pca_obs = mpca_scores(truth_real, :);
pca_scores_well = PredictionScoreWell(avail_real,:);


[rows, col] = size(pca_scores_well);
fprintf('The concatinated matrix has dimensions %d x %d.\n',rows, col);
  

end





