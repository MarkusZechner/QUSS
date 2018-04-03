function [ output_args ] = PlotDataCanonicalSpaceMoreD(D,DObs,Type,FontSize)
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
NumPlots=3;
MaxPlotPerRow=3;
NumDimensions=size(H,2);

for i=1:NumDimensions
    if mod(i,NumPlots)==1
        h=figure('Units', 'normalized', 'Position', [0,0,1,1]);
    end
    
    subplot(NumPlots/MaxPlotPerRow,MaxPlotPerRow,mod(i-1,NumPlots)+1);
    hold on;
    scatter(D(:,i),H(:,i),ScatterSize,'filled')
    plot([DObs(i),DObs(i)],[min(H(:,i)),max(H(:,i))],'r-',...
        'LineWidth',ObservedLineThickness);
    
    text(DObs(i) + abs(DObs(i))*0.25,min(H(:,i)) + ...
        abs(min(H(:,i)))*0.25,'d_{obs}','Fontweight','b','FontSize',FontSize);
    xlabel(['d_',num2str(i),'^' Type],'FontSize',FontSize);
    ylabel(['h_',num2str(i),'^' Type],'FontSize',FontSize);
    rho = corrcoef(D(:,i),H(:,i));
    title(['\rho = ' num2str(rho(1,2))],'FontSize',FontSize);
   
    set(gca,'FontSize',FontSize);
    axis square; axis tight;
    set(gcf,'color','w');
    %export_fig([FigureFolder 'ScreeMixedPCA' num2str(i)],'-m4','-transparent');
    
end

end



end

