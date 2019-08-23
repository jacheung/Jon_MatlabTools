
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses MATLAB's Breiman Random Forest Classifier to build decision trees.
% Input your design matrix X and the output class Y. 
%
%Will output predicted classes vs real classes for a random 30% of your
%data that wasnt trained on the trees. Though you can just use OOB
%prediction I kept it similar to how I was testing my logistic classifier.
%
% Will also output Mdl which are all the model parameters of the tree. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [predclasses, realclasses, motorClasses ,Mdl] = uber_baggingClassifier(DmatX,DmatY,motorX)

%biased sampling for training and testing for each class 
rando = randperm(length(DmatX));
tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
motorX = motorX(rando);
g1 = find(tmpDmatY == 1);
g2 = find(tmpDmatY == 2);
g1sidx = round(numel(g1)*.7);
g2sidx = round(numel(g2)*.7);
trainSet = [g2(1:g2sidx) ; g1(1:g1sidx)];
testSet = [1:length(DmatX)]';
testSet(trainSet)=[];

%Random Forests
numTrees = 50;
Mdl = TreeBagger(numTrees,tmpDmatX(trainSet,:),tmpDmatY(trainSet),'OOBPrediction','On',...
    'Method','classification','SampleWithReplacement','On','OOBPredictorImportance','On');

[fit,scores] =predict(Mdl,tmpDmatX(testSet,:));
predclasses= [str2double(fit) scores(:,1)];
realclasses = tmpDmatY(testSet);
motorClasses = motorX(testSet);

testerror = mean(predclasses(:,1)==realclasses);
disp(['Single Bagger Test Error = ' num2str(testerror*100)])

oobErrorBaggedEnsemble= oobError(Mdl);



%% plotting functions 
% figure (28);clf;subplot(2,2,[1 2]);plot(1-oobErrorBaggedEnsemble','k','linewidth',2)
% hold on; plot(mean(1-oobErrorBaggedEnsemble),'k','linewidth',2)
% xlabel 'Number of grown trees';
% ylabel 'OOB classification accuracy';
% 
% figure (28);subplot(2,2,3)
% bar(Mdl.OOBPermutedVarDeltaError)
% xlabel('Feature Index')
% ylabel('OOB Feature Importance')
% set(gca,'xlim',[0+.5 size(DmatX,2)+.5],'xticklabel',{'Theta','Counts','Pre 1 Lick','Pre 2 Lick','Pre 3 Lick'})
% xtickangle(45)
% 
% figure (28);subplot(2,2,4)
% gPosition = find(strcmp('1',Mdl.ClassNames));
% [Yfit,Sfit] = oobPredict(Mdl);
% [fpr,tpr,~,AUC] = perfcurve(Mdl.Y,Sfit(:,gPosition),'1');
% hold on;plot(fpr,tpr)
% xlabel('False Positive Rate')
% ylabel('True Positive Rate')
% text(.6,.5,['AUC = ' num2str(AUC)])


% plot(oobMeanMargin(Mdl));
% xlabel('Number of Grown Trees')
% ylabel('Out-of-Bag Mean Classification Margin')



% finbag = zeros(1,Mdl.NTrees);
% for t=1:Mdl.NTrees
%     finbag(t) = sum(all(~Mdl.OOBIndices(:,1:t),2));
% end
% finbag = finbag / size(DmatX,1);
% figure
% plot(finbag)
% xlabel('Number of Grown Trees')
% ylabel('Fraction of In-Bag Observations')





% [fpr,accu,thre] = perfcurve(Mdl.Y,Sfit(:,gPosition),'1','YCrit','Accu');
% figure(20)
% plot(thre,accu)
% xlabel('Threshold for ''good'' Returns')
% ylabel('Classification Accuracy')

%     view(Mdl.Trees{2},'Mode','graph')



%     figure(9);clf;plot(oobErrorBaggedEnsemble','color',[.8 .8 .8])
%     hold on; plot(mean(oobErrorBaggedEnsemble),'k','linewidth',2)
%     xlabel 'Number of grown trees';
%     ylabel 'Out-of-bag classification error';

