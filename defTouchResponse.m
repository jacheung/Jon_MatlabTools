%code used to define touch response. Touch response is the neural response
%post touch that exceeds mean + xCI 

confidenceThreshold = 0.99; 
rc = numSubplots(length(U));
figure(3000);clf

for rec=1:length(U)
    array = U{rec};
    window = [-25:50];
    touchIdx= [find(array.S_ctk(9,:,:)==1) ; find(array.S_ctk(12,:,:)==1)];
    spks = squeeze(array.R_ntk);
    
    blIdx = window(1:find(window==0));
    
    touchSpks = spks(touchIdx+window);
    touchSpksShuff = spks(touchIdx+blIdx);
    
    %calculating x% confidence interval
    SEM = nanstd(touchSpksShuff(:))./ sqrt(sum(~isnan(touchSpksShuff(:,1))));
    ts = tinv(confidenceThreshold,sum(~isnan(touchSpksShuff(:,1))));
    CI = SEM.*ts;
    
    touchResponse = smooth(mean(touchSpks));
    excitThreshold = mean(touchSpksShuff(:)) + CI;
    inhibThreshold = mean(touchSpksShuff(:)) - CI;
    
    %Excitatory threshold defined as the period > mean+x%CI
    excitthreshIdx = touchResponse'>excitThreshold;
    excitthreshIdx(1:find(window==0))=0;
    
    %Inhibitory threshold defined as the period > mean-x%CI
    inhibthreshIdx = touchResponse'<inhibThreshold;
    inhibthreshIdx(inhibthreshIdx==1)= -1;
    inhibthreshIdx(1:find(window==0))=0;

    figure(3000);subplot(rc(1),rc(2),rec)
    hold on; scatter(window,touchResponse*1000,'k')
    hold on; plot(window,ones(length(window),1).* excitThreshold .* 1000,'r-.')
    
    %Defining touch responses as a period between two points that are less
    %than 5ms apart. 
    if sum(excitthreshIdx)>0
        tps = window(excitthreshIdx);
        if ~isempty(tps)
            startPoint = tps(1);
            endPoint = tps(find(diff(diff(tps)<5)==-1,1,'first'));
            if isempty(endPoint)
                endPoint = tps(end);
            end
            U{rec}.meta.responseWindow=[startPoint endPoint];
            hold on; scatter(window(startPoint+find(window==0):endPoint+find(window==0)),touchResponse(startPoint+find(window==0):endPoint+find(window==0))*1000,'b','filled')
        end
    end

end