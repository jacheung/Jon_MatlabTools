% Function returns indices of each class for stratified training and testing
%
% Input vector of prediction classes and proportion for training
% Output cell for each prediction class

function [trainIdx,testIdx] = stratifiedIndices(DmatY,proportionTrain)

if nargin<2
    error('must input proportion to set aside for train (e.g. .7)')
end

classes = unique(DmatY);

trainIdx = cell(1,numel(classes)); 
testIdx = cell(1,numel(classes)); 

for d = 1:numel(classes)
    classidx = find(DmatY==classes(d));
    shuffIdx = classidx(randperm(length(classidx)));
    
    trainIdx{d} = shuffIdx(1:round(length(classidx)*proportionTrain));
    testIdx{d} = setdiff(shuffIdx,trainIdx{d});
end


