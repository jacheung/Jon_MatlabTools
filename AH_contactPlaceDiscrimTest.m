numCell = length(L);

PLickTouchGo = zeros(numCell,1);
PLickNoTouchGo = zeros(numCell,1);
PLickTouchNoGo = zeros(numCell,1);
PLickNoTouchNoGo = zeros(numCell,1);

numLickTouchGo = zeros(numCell,1);
numLickNoTouchGo = zeros(numCell,1);
numLickTouchNoGo = zeros(numCell,1);
numLickNoTouchNoGo = zeros(numCell,1);

numTrialHit = zeros(numCell,1); % Lick Go
numTrialMiss = zeros(numCell,1);
numTrialCR = zeros(numCell,1);
numTrialFA = zeros(numCell,1); % Lick Nogo
numTrialGo = zeros(numCell,1); % total go trials    
numTrialNoGo = zeros(numCell,1); % total no-go trials

for i = 1 : numCell
    % refresh temp matrices
    tempTouchGo = [];
    tempNoTouchGo = [];
    tempTouchNoGo = [];
    tempNoTouchNoGo = [];
    
    % loading data
    if i < 10
        cellname = strcat('Cell0',num2str(i),'*.mat');
    else
        cellname = strcat('Cell',num2str(i),'*.mat');
    end
    cellname2 = dir(cellname);
    load(cellname2.name);
    
    % find out time series length of each trial
    lengthEachTrial = max(find(c.timeSeriesArrayHash.value{1}.trial == min(c.timeSeriesArrayHash.value{1}.trial)));
    
    % find out trial indices for go and no-go trials
    trialHit = find(c.trialTypeMat(1,:) == 1); % Lick Go
    trialMiss = find(c.trialTypeMat(2,:) == 1);
    trialCR = find(c.trialTypeMat(3,:) == 1);
    trialFA = find(c.trialTypeMat(4,:) == 1); % Lick Nogo
    trialGo = sort([trialHit, trialMiss]); % total go trials    
    trialNoGo = sort([trialCR, trialFA]); % total no-go trials
    
    numTrialHit(i) = length(trialHit); % Lick Go
    numTrialMiss(i) = length(trialMiss);
    numTrialCR(i) = length(trialCR);
    numTrialFA(i) = length(trialFA); % Lick Nogo
    numTrialGo(i) = length(trialGo); % total go trials    
    numTrialNoGo(i) = length(trialNoGo); % total no-go trials

    % find out trial indices for touched during go, not touched during go,
    % touched during no-go, and not touched during no-go
    for j = trialGo
        touchInd = sum(c.timeSeriesArrayHash.value{1,1}.valueMatrix(7,(j-1)*lengthEachTrial+1:j*lengthEachTrial)); % touch_onset, if 0, there was no touch
        if touchInd == 0 % not touched
            tempNoTouchGo = [tempNoTouchGo, j];
        else
            tempTouchGo = [tempTouchGo,j];
        end
    end
    for j = trialNoGo
        touchInd = sum(c.timeSeriesArrayHash.value{1,1}.valueMatrix(7,(j-1)*lengthEachTrial+1:j*lengthEachTrial)); % touch_onset, if 0, there was no touch
        if touchInd == 0
            tempNoTouchNoGo = [tempNoTouchNoGo, j];
        else
            tempTouchNoGo = [tempTouchNoGo, j];
        end
    end
    
    % numbers for calculating the probability
    numLickTouchGo(i) = sum(ismember(trialHit, tempTouchGo));
    numLickNoTouchGo(i) = sum(ismember(trialHit, tempNoTouchGo));
    numLickTouchNoGo(i) = sum(ismember(trialFA, tempTouchNoGo));
    numLickNoTouchNoGo(i) = sum(ismember(trialFA, tempNoTouchNoGo));
    
    PLickTouchGo(i) = numLickTouchGo(i)/length(tempTouchGo);
    PLickNoTouchGo(i) = numLickNoTouchGo(i)/length(tempNoTouchGo);
    PLickTouchNoGo(i) = numLickTouchNoGo(i)/length(tempTouchNoGo);
    PLickNoTouchNoGo(i) = numLickNoTouchNoGo(i)/length(tempNoTouchNoGo);

end

figure, hold on, scatter(zeros(numCell,1)+1, PLickTouchGo), scatter(zeros(numCell,1)+2, PLickNoTouchGo), 
scatter(zeros(numCell,1)+3, PLickTouchNoGo), scatter(zeros(numCell,1)+4, PLickNoTouchNoGo)
save('contactPlaceDiscrim.mat', 'P*', 'numLick*', 'trial*')
