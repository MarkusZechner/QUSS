% function [h_reconstructed] = UndoCanonicalAndPCA(h_c_post,B, h_f,Time,...
% PriorPredictionFPCA)
%
% UndoCanonicalFunctional Undo the canonical and functional data analysis
% transformations
%   Posterior samples are generated in the canonical space, we need to undo
%   the canonical transform as well as the functional data transformations
%   to return the samples back into the time domain.
%
% Inputs:
%   h_c_post: (NReals x NCanonicalDim) posterior samples
%   B: (NOriginalDim x NCanonicalDim) Transformation matrix from CCA
%   Time: (NTimeSteps) Vector of time steps that we need to project back onto
%   PriorPredictionFPCA: FPCA struct that contains the eigenfunctions need
%   to undo that transformation
%
% Outputs:
%   h_reconstructed: (NReals X NOriginalDim) Reconstructed posterior
%   realizations

function [h_reconstructed] = UndoCanonicalAndPCA(h_c_post,B, h_f,...
    PriorPredictionPCA, nEigenvalues)

% Get number of posterior samples
NumPosteriorSamples = length(h_c_post);

% Undo the canonical correlation analysis
% was named HpostCoef before but not accurate I think??
HpostScore = h_c_post'*pinv(B)+repmat(mean(h_f,1)',...
    1,NumPosteriorSamples)';


h_reconstructed = HpostScore * PriorPredictionPCA{1}.coeff(:,1:nEigenvalues)';
h_reconstructed = bsxfun(@plus, h_reconstructed, PriorPredictionPCA{1}.mu);



end