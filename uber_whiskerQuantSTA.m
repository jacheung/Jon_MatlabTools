%% Whisking quantification
lookbackwndw=50;

for d = 1:length(U)
    masks =  assist_touchmasks(U{d});
    responses =  squeeze(U{d}.R_ntk);
    responses=responses.*masks.touch;
    responses(1:lookbackwndw,:)=0;
    
    responseIdx = find(responses==1);
    figure(380);clf
    stimVals = [1 2 3 4 5 6];
    for g = 1:length(stimVals)
        
        stimuli = squeeze(U{d}.S_ctk(stimVals(g),:,:));
        
        raw = zeros(length(responseIdx),lookbackwndw);
        avgs = zeros(1,lookbackwndw);
        for n = 1:length(responseIdx)
            stimIdx = responseIdx(n)-lookbackwndw+1 : responseIdx(n);
            
            if sum(isnan(stimuli(stimIdx)))>0
                disp('Skipping spike b/c NaN value')
            else
                raw(n,:) = stimuli(stimIdx);
                avgs=avgs+stimuli(stimIdx);
            end
        end
        
        
        figure(380);
        subplot(2,3,g)
        boundedline(1:lookbackwndw,mean(raw),std(raw));
        set(gca,'xtick',0:10:50,'xticklabel',-50:10:0,'xlim',[0 50])
    end
    pause
end

