function [predclasses, realclasses, motorClasses ,Mdl1, Mdl2] = uber_DoublebaggingClassifier(DmatX,DmatY,motorX,uber,splitmethod)

%biased sampling for training and testing for each class 
notouchidx= find(DmatX(:,2) == 0);
Uranges = [uber.meta.ranges];
splitmotors = [mean(Uranges) - 15000 mean(Uranges) + 15000];
motordefidx = intersect(find(motorX>=splitmotors(1)),find(motorX<=splitmotors(2)));

if strcmp(splitmethod,'notouches')
    splitidx = notouchidx; %this will define what type of splitting b/t the two baggers
elseif strcmp(splitmethod,'motordefined')
    splitidx = motordefidx;
else 
    error('Choose split method: notouches or motordefined')
end

split1X = DmatX(splitidx,:);
split1Y = DmatY(splitidx,:);
split1Motor = motorX(splitidx,:);

DmatX(splitidx,:) =[];
DmatY(splitidx,:) =[];
motorX(splitidx,:) = [];


%% Building Random Forest for all trials WITH touches 
rando = randperm(size(DmatX,1));
tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
motorX = motorX(rando);
g1 = find(tmpDmatY == 1);
g2 = find(tmpDmatY == 2);
g1sidx = round(numel(g1)*.7);
g2sidx = round(numel(g2)*.7);
trainSet = [g2(1:g2sidx) ; g1(1:g1sidx)];
testSet = [1:size(DmatX,1)]';
testSet(trainSet)=[];

%Random Forests
numTrees = 50;
Mdl1 = TreeBagger(numTrees,tmpDmatX(trainSet,:),tmpDmatY(trainSet),'OOBPrediction','On',...
    'Method','classification','SampleWithReplacement','On','OOBPredictorImportance','On');

[fit,scores] =predict(Mdl1,tmpDmatX(testSet,:));
predclasses = [str2double(fit) scores(:,1)];
realclasses = tmpDmatY(testSet);
motorClasses = motorX(testSet);

oobErrorBaggedEnsemble= oobError(Mdl1);
figure (28);clf;subplot(2,2,[1 2]);plot(1-oobErrorBaggedEnsemble','k','linewidth',2)
ylabel 'OOB classification accuracy';


%% Building Random Forest for all trials WITHOUT touches 
rando = randperm(size(split1X,1));
mtDmatX=split1X(rando,:);mtDmatY=split1Y(rando,:);
mtmotorX = split1Motor(rando);
g1 = find(mtDmatY == 1);
g2 = find(mtDmatY == 2);
g1sidx = round(numel(g1)*.7);
g2sidx = round(numel(g2)*.7);
trainSet = [g2(1:g2sidx) ; g1(1:g1sidx)];
testSet = [1:size(split1X,1)]';
testSet(trainSet)=[];

%Random Forests
numTrees = 50;
Mdl2 = TreeBagger(numTrees,mtDmatX(trainSet,:),mtDmatY(trainSet),'OOBPrediction','On',...
    'Method','classification','SampleWithReplacement','On','OOBPredictorImportance','On');

[fit,scores] =predict(Mdl2,mtDmatX(testSet,:));
predclasses= [predclasses ; str2double(fit) scores(:,1)];
realclasses = [realclasses ; mtDmatY(testSet)];
motorClasses = [motorClasses; mtmotorX(testSet)];

testerror = mean(predclasses(:,1)==realclasses);
disp(['Split Bagger Total Test Error = ' num2str(testerror*100)])

oobErrorBaggedEnsemble= oobError(Mdl2);
figure (28);hold on;subplot(2,2,[3 4]);plot(1-oobErrorBaggedEnsemble','k','linewidth',2)
xlabel 'Number of grown trees';
ylabel ('OOB classification accuracy');