function [F1scores] = F1scorev2(predictions,actual,numClasses)
%% OUTPUTS PREDICTIONS which is is 3 columns (goF1, nogo F1)



F1scores = zeros(1,numClasses);

for i = 1:numClasses

PredClass = (predictions==i);
ActualClass = (actual==i); %all true positives 

TP = sum(PredClass(ActualClass==1)==1);
FP = sum(PredClass((ActualClass==0))==1);
FN = sum(PredClass(ActualClass== 1)==0);

P = TP./(TP+FP);
R = TP./(TP+FN);

% 
% TP = sum(predictions(ActualPos)==i);
% FN = sum(~(predictions(ActualPos)==i));
% predPos = sum(predictions==i);

% P = TP./predPos; %precision defined as TRUEPOS/predicted pos
% 
% R = TP./(TP+FN); %recall defined as TRUEPOS/ACTUALPOS
    
F1s = 2 .* (P.*R)./(P+R);

F1scores(i) = F1s;
end 


%GROUP F1
TP = sum(predictions(actual==1)==1);
FP = sum(predictions((actual==2))==1);
FN = sum(predictions(actual == 1)==2);


P = TP./(TP+FP);
R = TP./(TP+FN);
F1all = 2 .* (P.*R)./(P+R);
F1scores = [F1scores F1all];