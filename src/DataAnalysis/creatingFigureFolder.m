function [ FigureFolder ] = creatingFigureFolder(SaveFig, BasePath, Case, Trial, WellName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if SaveFig == 'On'
    
    if (nargin > 4)
        Name = WellName{1};
        FigureFolder = [BasePath '\' Case '\' Trial '_' Name '\'];
    else
        FigureFolder = [BasePath '\' Case '\' Trial '\'];
    end
    
    
    % check if it is existing, if not make a new one
    if exist(FigureFolder,'dir') ~= 7
        mkdir(FigureFolder);
    end
    
end

