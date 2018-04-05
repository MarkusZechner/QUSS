% This script uses PCA components of data and forcast and 
% identifies the sensitivity between them.
% Well Data are the parameters here
% the PCA components of the forecast are the objective function


% could check the influence on the all the well on the field.
% shouldnt we see the largest producers being most important?

close all; clear all; clc;
addpath('../src/DataAnalysis/');
addpath('../src/thirdparty/export_fig/');
addpath('../../DGSA/');

prior_path = '../data/evidential/PriorData.mat';
load(prior_path);
prior_path = '../data/evidential/PriorPrediction.mat';
load(prior_path);

% Saving Figures
BasePath = 'C:\Users\Markus Zechner\Documents\GitHub\QUSS\figures\DataAnalysis';
Case = 'OMV8TH_largerPrior';
% Trial is a subfolder of Case
Trial = 'Test2';

SaveFig = 'On';
[ FigureFolder ] = creatingFigureFolder(SaveFig, BasePath, Case, Trial);

%% Performing PCA on the time signal directly (no B-splines)on both d and h
% the EigenTolerance in the next 2 functions is only used to demonstrate 
% the reconstruction of the signal but is not used to remove the ones 
%  describing less varaince, this is done when the wells are concatinated
rmpath('../src/thirdparty/fda_matlab/');
PlotLevel = 1; % 1=yes, 2=no

MinEigenValues = 3; EigenTolerance = 0.97;

PriorDataPCA = ComputePCADataAnalysis(PriorData, EigenTolerance, MinEigenValues, PlotLevel, ...
    'Data', FigureFolder );

PriorPredictionPCA = ComputePCADataAnalysis(PriorPrediction, EigenTolerance, MinEigenValues, PlotLevel, ...
    'Prediction', FigureFolder );

% Concatinating the matrix

NumberOfEigenValuesData = 2;
NumberOfEigenValuesPrediction = 31;
TruthRealization = 1;
Normalize = true;
WellNumber = 68;
FontSize = 20;

% concatinate PCA scores of all wells and remove truth
[pca_scores_data] = ConcatinatePCAMatrix(PriorDataPCA,TruthRealization,...
    NumberOfEigenValuesData, Normalize);

% if more than 1 PCA Component is used, the wellnames have to be adjusted

[ NewWellNames ] = MakeNewWellNames( PriorData, NumberOfEigenValuesData );
PriorData.ObjNames = NewWellNames;
% ploting the Data of the selected well
PlotWellResponses(PriorData,TruthRealization,WellNumber, FontSize)

% ploting the Prediction of the selected well
PlotWellResponses(PriorPrediction,TruthRealization,WellNumber, FontSize)

% get PCA scores for a specific well and remove the truth
[pca_scores_well] = getPCAscoresPrediction(PriorPredictionPCA, ...
    TruthRealization, NumberOfEigenValuesPrediction,WellNumber);

% computing the distance between the prediction responses (selected well)

Distance = squareform(pdist(pca_scores_well));

DGSA = {};
DGSA.ParametersNames = PriorData.ObjNames;
DGSA.ParametersValues = pca_scores_data;
DGSA.D = Distance;
DGSA.N = size(PriorData.data,1)-1;

save('C:\Users\Markus Zechner\Documents\GitHub\QUSS\data\evidential\DistanceMAllComp2.mat','DGSA');

