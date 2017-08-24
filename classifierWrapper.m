% Function for wrapping necessary variables for the logistic classifier.
% Takes the uber array and a couple other functions to output a struct with
% all desired predecision touch variables 

%INPUT: uber array
%OUTPUT: V array which is a struct for all predecision touch variables

function [V] = classifierWrapper(uberarray)

for rec = 1:length(uberarray)
    
    array=uberarray{rec};
    
    V(rec).varNames = {'theta','velocity','amplitude','setpoint','phase','deltaKappa'};
    
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    varx=[1:6];
    for f=1:length(varx)
        [~, ~,V(rec).var.hit{f}, V(rec).var.miss{f}, V(rec).var.FA{f}, V(rec).var.CR{f},~,~] = assist_vardistribution(array,varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
        if f == 1
            V(rec).touchNum.hit=cellfun(@numel,V(rec).var.hit{1}); % count number of predecision touches using theta vals
            V(rec).touchNum.miss=cellfun(@numel,V(rec).var.miss{1});
            V(rec).touchNum.FA=cellfun(@numel,V(rec).var.FA{1});
            V(rec).touchNum.CR=cellfun(@numel,V(rec).var.CR{1});
        end
        V(rec).var.hit{f}   =      cell2mat(V(rec).var.hit{f}); %predecision variables for hit trials
        V(rec).var.miss{f}  =      cell2mat(V(rec).var.miss{f});
        V(rec).var.FA{f}    =      cell2mat(V(rec).var.FA{f});
        V(rec).var.CR{f}    =      cell2mat(V(rec).var.CR{f});     
    end
    
    V(rec).trialNums.matNames = {'hit','miss','FA','CR','licks'};
    V(rec).trialNums.matrix = zeros(5,array.k);
    V(rec).trialNums.matrix(1,find(array.meta.trialType == 1 & array.meta.trialCorrect ==1))= 1; %hit trials
    V(rec).trialNums.matrix(2,find(array.meta.trialType == 1 & array.meta.trialCorrect ==0))= 1; %miss trials
    V(rec).trialNums.matrix(3,find(array.meta.trialType == 0 & array.meta.trialCorrect ==0))= 1; % FA trials
    V(rec).trialNums.matrix(4,find(array.meta.trialType == 0 & array.meta.trialCorrect ==1))= 1; %CR trials
    V(rec).trialNums.matrix(5,:) = sum(V(rec).trialNums.matrix([1 3],:));

    %figuring out previous TT and whether it was lick or not
    lix{1} = [0 V(rec).trialNums.matrix(5,:)]; %padded with 0 to account for indexing
    lix{2} = [0 0 V(rec).trialNums.matrix(5,:)];
    lix{3} = [0 0 0 V(rec).trialNums.matrix(5,:)];
        
    hx_prevT = find(V(rec).trialNums.matrix(1,:)==1); % using current T num b/c licks shifted by padded 0
    mx_prevT = find(V(rec).trialNums.matrix(2,:)==1);
    FAx_prevT = find(V(rec).trialNums.matrix(3,:)==1); % using current T num b/c licks shifted by padded 0
    CRx_prevT = find(V(rec).trialNums.matrix(4,:)==1); % using current T num b/c licks shifted by padded 0
    
    V(rec).licks.oneT.hit =lix{1}(hx_prevT);
    V(rec).licks.oneT.miss =lix{1}(mx_prevT);
    V(rec).licks.oneT.FA =lix{1}(FAx_prevT);
    V(rec).licks.oneT.CR =lix{1}(CRx_prevT);
    
    V(rec).licks.twoT.hit =lix{2}(hx_prevT);
    V(rec).licks.twoT.miss =lix{2}(mx_prevT);
    V(rec).licks.twoT.FA =lix{2}(FAx_prevT);
    V(rec).licks.twoT.CR =lix{2}(CRx_prevT);
    
    V(rec).licks.threeT.hit =lix{3}(hx_prevT);
    V(rec).licks.threeT.miss =lix{3}(mx_prevT);
    V(rec).licks.threeT.FA =lix{3}(FAx_prevT);
    V(rec).licks.threeT.CR =lix{3}(CRx_prevT);
    
end