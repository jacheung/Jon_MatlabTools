function [F1scores] = F1score(predictions,actual,numClasses)

F1scores = zeros(1,numClasses);

for i = 1:numClasses

ActualPos =(actual==i); %all true positives 

TP = sum(predictions(ActualPos)==i);
FN = sum(~(predictions(ActualPos)==i));
predPos = sum(predictions==i);

P = TP./predPos; %precision defined as TRUEPOS/predicted pos

R = TP./(TP+FN); %recall defined as TRUEPOS/ACTUALPOS
    
F1s = 2 .* (P.*R)./(P+R);

F1scores(i) = F1s;
end 