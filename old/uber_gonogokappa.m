figure(2800);clf
for rec = 1:length(BV)
    clearvars -except BV rec
    array = BV{rec};
    gos = find(array.meta.trialType==1);
    nogos = find(array.meta.trialType==0);
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    
    groupKappgo = [];
    for d = 1:length(prelixGo)
        if ~isempty(prelixGo{d})
            phase = array.S_ctk(5,:,gos(d));
            kapps = array.S_ctk(6,:,gos(d));
            valIdx = repmat(prelixGo{d}',1,101) + repmat(0:100,length(prelixGo{d}),1);
%             tmp = phase(valIdx);
%             tmp2 = find(tmp(:,1)<0);
%             if ~isempty(tmp2)   
                groupKappgo = [groupKappgo ; kapps(valIdx)];
%             end
        end
    end
    
    groupKappnogo = [];
    for d = 1:length(prelixNoGo)
        if ~isempty(prelixNoGo{d})
            kapps = array.S_ctk(6,:,nogos(d));
            valIdx = repmat(prelixNoGo{d}',1,101) + repmat(0:100,length(prelixNoGo{d}),1);
            groupKappnogo = [groupKappnogo ; kapps(valIdx)];
        end
    end
    
    x = randi(50, 1, 100);                      % Create Data
SEM = nanstd(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.05  0.95],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals
    
    % 95% CI
    for i = 1:size(groupKappgo,2)
        x=groupKappgo(:,i);
        SEM = nanstd(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.05  0.95],length(x)-1);      % T-Score
        Tcibingo([2 1],i) = nanmean(x) + ts*SEM;   %confidence intervals
        
        z=groupKappnogo(:,i);
        SEM = nanstd(z)/sqrt(length(z));               % Standard Error
        ts = tinv([0.05  0.95],length(z)-1);      % T-Score
        Tcibinnogo([2 1],i) = nanmean(z) + ts*SEM;   
    end
    
    figure(2800);subplot(2,5,rec)
    plot(0:100,nanmean(groupKappgo),'b');
    hold on; plot(0:100,Tcibingo(1,:),'-.b')
    hold on; plot(0:100,Tcibingo(2,:),'-.b')
    plot(0:100,nanmean(groupKappnogo),'r');
    hold on; plot(0:100,Tcibinnogo(1,:),'-.r')
    hold on; plot(0:100,Tcibinnogo(2,:),'-.r')
end

%
% figure;boundedline(0:100,mean(groupKappgo),std(groupKappgo),'b')
% hold on;boundedline(0:100,mean(groupKappnogo),std(groupKappnogo),'r')

