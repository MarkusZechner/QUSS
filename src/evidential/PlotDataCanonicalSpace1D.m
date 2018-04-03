function [ output_args ] = PlotDataCanonicalSpace1D(D,DObs,Type,FontSize, FigureFolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Inputs:
%   D: (NReals x NDDim) Data variable
%   DObs: (1 x NDim) Observed data
%   Type: c or f for canonical space or functional space
%   FontSize: (int) Font size for display [optional]
% Return:
%   h: handle to figure
%
% Author: Markus Zechner (mzechner@stanford.edu)
% Date: May 2nd 2016

if (nargin < 4)
    FontSize=12;
end

ScatterSize=50;
ObservedLineThickness=3;


i=1; % only plotting the firs dimension

figure;
h1=scatter(D(:,i),D(:,i+1),ScatterSize,[0.5 0.5 0.5]);
hold on
h2=scatter(DObs(i),DObs(i+1),65,'r','filled');
legend([h1(1),h2(1)],'Prior','d_{obs}','Location','northwest');

xlabel(['d_',num2str(i),'^' Type],'FontSize',FontSize);
ylabel(['d_',num2str(i+1),'^' Type],'FontSize',FontSize);

set(gca,'FontSize',FontSize);
axis square; axis tight;
set(gcf,'color','w');
export_fig([FigureFolder 'DataCanonicalSpace2D'],'-m4','-transparent');



end


