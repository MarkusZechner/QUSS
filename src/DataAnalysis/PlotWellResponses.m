function [ h ] = PlotWellResponses(DataStruct,TruthRealization,WellNumber, FontSize)
%PLOTRESPONSES: Plots input response for data and prediction variables
%
% Inputs:
%   HistoricalStruct: Struct containing historical data
%   ForecastStruct: Struct containing forecast data
%   TruthRealization: Index of realization taken to be the truth
%   FontSize: Font size in plots
%
% Author: Lewis Li (lewisli@stanford.edu)
% Original Date: March 4th 2016
% Last Updated: September 26th 2016

if (nargin < 3)
    FontSize=12;
end

NumRealizations = size(DataStruct.data,1);

figure;


for j = 1:NumRealizations
    h1 = plot(DataStruct.time, DataStruct.data(j,:,WellNumber)',...
        'color',[0.5 0.5 0.5]);
    hold on
end
xlabel('Time (days)','FontSize',FontSize);
ylabel(DataStruct.name,'FontSize',FontSize);


if (isstruct(TruthRealization))
    h2 = plot(TruthRealization.time, ...
        TruthRealization.data(:,:,WellNumber)','r','LineWidth',3);
    hlegend = legend([h1,h2],'Prior','Observed');
    set(hlegend,'Location','southwest');
elseif TruthRealization>0
    h2 = plot(DataStruct.time, ...
        DataStruct.data(TruthRealization,:,WellNumber)','r','LineWidth',3);
    hlegend = legend([h1,h2],'Prior','Observed');
    set(hlegend,'Location','southwest');
end

title([ DataStruct.type ': ' DataStruct.ObjNames{WellNumber}],...
    'FontSize',FontSize); axis square;

axis tight;
set(gcf,'color','w');
set(gca,'FontSize',FontSize);

end



