function [ NewNames ] = MakeNewWellNames( PriorData, NumberOfEigenValuesData  )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here


NewNames = {};

Names = PriorData.ObjNames;

NbNames = size(Names, 2);

k=1;
for i=1:NbNames
    
    for j=1:NumberOfEigenValuesData
        
        NewNames{k} = [Names{i} '_C_' num2str(j)];
        
        k = k + 1;
        
    end
end





end

