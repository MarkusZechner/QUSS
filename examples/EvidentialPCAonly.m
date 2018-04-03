
% This script is an example for:
% Wells - Field or 
% Field - Field
% using PCA as the dim reduction method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Evidential.m
% Tutorial file on how to use Evidential Learning for forecasting future
% prediction rates from historical production rates.
%
% Author: Lewis Li (lewisli@stanford.edu)
% Original Date: December 30th 2016
% Last Modified: Janurary 30th 2017

close all; clear all; clc;
addpath('../src/evidential/');
addpath('../src/thirdparty/export_fig/');

prior_path = '../data/evidential/PriorData.mat';
load(prior_path);
prior_path = '../data/evidential/PriorPrediction.mat';
load(prior_path);

% Saving the images
BasePath = 'C:\Users\Markus Zechner\Documents\GitHub\QUSS\figures\EvidentialLearning';
Case = 'OMV8TH_largerPrior';
% Trial is a subfolder of Case
Trial = 'FieldToField_PCAonly_EigenTol097';

% creating the figure folder
SaveFig = 'On';
[ FigureFolder ] = creatingFigureFolder(SaveFig, BasePath, Case, Trial);



%% creating a field response

% Data
PriorRates = PriorData.data;
PriorRatesObserved = PriorRates(1,:,:);
PriorRatesSimulated = PriorRates(2:end,:,:);

PriorFieldObserved = sum(PriorRatesObserved,3);
PriorFieldSimulated = sum(PriorRatesSimulated,3);

DataNewPrior = [PriorFieldObserved;PriorFieldSimulated];
PriorData.data = DataNewPrior;
PriorData.ObjNames = {'Field'};

% Forecast
PredictionRates = PriorPrediction.data;
PredictionRatesObserved = PredictionRates(1,:,:);
PredictionRatesSimulated = PredictionRates(2:end,:,:);

PredictionFieldObserved = sum(PredictionRatesObserved,3);
PredictionFieldSimulated = sum(PredictionRatesSimulated,3);

DataNewPrediction = [PredictionFieldObserved;PredictionFieldSimulated];
PriorPrediction.data = DataNewPrediction;
PriorPrediction.ObjNames = {'Field'};

%% Plot input responses
FontSize = 20;


% Set aside a realization to use as the "truth"
TruthRealization = 1; 
NumPriorRealizations=length(PriorData.data);
AvailableRealizations = setdiff(1:NumPriorRealizations,TruthRealization);

%PlotResponses(PriorData,TruthRealization,FontSize);
%PlotPriorResponses(PriorPrediction,TruthRealization,FontSize);

%% Dimension Reduction On Both Data and Prediction Variables

% Minimum number of eigenvalues and of 'explained variance' to keep after dim
% reduction
% IMPORTANT: this is only used to demonstrate how well the dimension
% reduction worked. ALL scores are transferd to mixedPCA. In the 
% mixed PCA then only the highest eigenvalues are kept!
MinEigenValues = 3; EigenTolerance = 0.97;

% Performing PCA on the time signal directly (no B-splines)on both d and h
rmpath('../src/thirdparty/fda_matlab/');
PlotLevel = 1; % 1=yes, 2=no
PriorDataPCA = ComputePCA(PriorData, EigenTolerance, MinEigenValues, PlotLevel, ...
    'Data', FigureFolder );

PriorPredictionPCA = ComputePCA(PriorPrediction, EigenTolerance, MinEigenValues, PlotLevel, ...
    'Prediction', FigureFolder );

% Perform Mixed PCA on FPCA components for d

[d_f, dobs_fpca] = MixedPCAnoBsplines(PriorDataPCA, TruthRealization, EigenTolerance, FigureFolder);

% Get number of FPCA components to keep for h
nEigenvalues = GetNumEigenvalues(PriorPredictionPCA{1}, MinEigenValues, ...
    EigenTolerance);
h_f = PriorPredictionPCA{1}.score(AvailableRealizations,1:nEigenvalues);

% Plot prior models in functional space
PlotLowDimModels(d_f,h_f,dobs_fpca,'f',FontSize);

%% Apply CCA transformation and compute posterior distribution in canonical space

% Compute canonical transformation on the functional components
[A, B, ~, d_c,h_c] = canoncorr(d_f,h_f);
dobs_c=(dobs_fpca-mean(d_f))*A;

% Plot prior models in canonical space
PlotLowDimModels(d_c,h_c,dobs_c,'c',FontSize);

% Plot the data in cannonical Space  including DObs
% this is to check whether they are gaussian
PlotDataCanonicalSpace1D(d_c,dobs_c,'c',FontSize, FigureFolder);

% Apply a normal score transform to the h_c
h_c_gauss = NormalScoreTransform(h_c,0);

% Find best linear fit between Dc and h_c_gauss
G = d_c'/h_c_gauss';

% Compute misfit covariance
DDiff= d_c'-G*h_c_gauss';
C_T = DDiff*DDiff'/length(d_c);

% Compute posterior using Linear Gaussian Regression Equations
C_H = cov(h_c_gauss);
mu_prior = mean(h_c_gauss)';
 mu_posterior = mu_prior + C_H*G'*pinv(G*C_H*G' + C_T)*(dobs_c'-G*mu_prior);
C_posterior = inv(G'*pinv(C_T)*G + inv(C_H));

%% Generate samples from the posterior in canonical space
%addpath('../src/thirdparty/fda_matlab/');
NumPosteriorSamples = 100;
h_c_post = SampleCanonicalPosterior(mu_posterior,C_posterior,...
    NumPosteriorSamples,h_c, FigureFolder );

% Undo the CCA and PCA transformations
h_reconstructed = UndoCanonicalAndPCA(h_c_post, B, h_f,...
    PriorPredictionPCA, nEigenvalues);

%% Compute and plot posterior responses and quantiles
% had to take away the ' is that OK??
[PriorQuantiles, PosteriorQuantiles] = ComputeQuantiles(...
    PriorPrediction.data, h_reconstructed);
PlotPosteriorSamplesAndQuantiles(PriorPrediction,TruthRealization, ...
    h_reconstructed,PriorQuantiles,PosteriorQuantiles, FigureFolder);


