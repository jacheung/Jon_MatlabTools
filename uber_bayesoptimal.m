
%% Build for a ghetto decoder.
%P(spks|stimuli) given as just a poisson distribution with lambda of the
%mean of all spks generated within the window for that angle bin
%P(stimuli|spks) given as just a normal distribution with mean of all
%stimuli that gave x number of spks.
%Decoder samples spikes for each bin based on P(spks|stimuli) and finds the
%P(bin|spks) using the normal distributions.




window = [8 25];
countThresh = 10;
trials2sample = 100;
neuronsample = 200;

figure(30);clf
figure(40);clf;
figure(50);clf;
figure(60);clf;
figure(80);clf;
figure(90);clf;
figure(100);clf;
clear decodeAcc
clear pstimgspks
clear resol
clear resolstd

for rec = [1 3 4 5 9 10 11 12 13]
    
    
    
    
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
    V.bounds = [{[-20:5:80]} {[-9750:125:9750]} {[-99.5:5:99.5]} {[-99.5:5:99.5]} {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
    
    
    for d = 1
        binrangetmp = V.bounds{d};
        binspace = (binrangetmp(end)-binrangetmp(1))/(numel(binrangetmp)-1);
        binrange = binrangetmp(1)+binspace/2:binspace:binrangetmp(end)-binspace/2;
        bloop = parula(size(binrange,2));
        
        %Build Variables
        var = squeeze(U{rec}.S_ctk(d,:,:));
        varAtTouchraw = var(touchIdx);
        varAtTouchraw = varAtTouchraw(~isnan(varAtTouchraw));
        touchIdx = touchIdx(~isnan(varAtTouchraw));
        spikesAtTouchraw = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(varAtTouchraw),1)),2);
       
        
        %remove any indices with only one example of spikes since we cant
        %do leave one out to test classifier with that
        if numel(find(spikesAtTouchraw == max(spikesAtTouchraw)))<=1
            spikesAtTouchraw(spikesAtTouchraw == max(spikesAtTouchraw))=[];
            varAtTouchraw(spikesAtTouchraw == max(spikesAtTouchraw))=[];
        end
        
        dropraw = randperm(length(varAtTouchraw));
        
        bayesclass = nan(length(dropraw),length(binrange)+2);
        for p = 1:length(dropraw)
            drop=dropraw(p);
            loo = [spikesAtTouchraw(drop) varAtTouchraw(drop)];
            spikesAtTouch = spikesAtTouchraw;
            varAtTouch = varAtTouchraw;
            spikesAtTouch(drop) = [];
            varAtTouch(drop)=[];

            [sorted, sortedBy ,binBounds]=binslin(varAtTouch,spikesAtTouch,'equalE',numel(V.bounds{d}),V.bounds{d}(1),V.bounds{d}(end));
            
            trims=[binrange' cell2mat(cellfun(@mean,sorted,'Uniformoutput',0)) cell2mat(cellfun(@std,sorted,'Uniformoutput',0))];
            counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
            mintouches=find(counttmp(:,1)>countThresh);
            
            %% building P(spks) using kernel
            
            spksx = [0:1:max(spikesAtTouch)];
            %pdfspks = normpdf(spksx,mean(spikesAtTouch),std(spikesAtTouch));
            h=histogram(spikesAtTouch,'normalization','probability');
            pdfspks = h.Values;
%             ker_vals = fitdist(spikesAtTouch,'kernel');
%             pdfspks = pdf(ker_vals,spksx);
            
            %% building P(stimuli) using kernel
            kervar_vals = fitdist(varAtTouch,'kernel');
            tmppdfstim = pdf(kervar_vals,binBounds(1):1:binBounds(end));
            
            bins = 1:binspace:binBounds(end)-binBounds(1);
            pdfstim = nan(1,length(bins));
            for m=1:length(bins)
                pdfstim(m) = sum(tmppdfstim(bins(m):bins(m)+binspace-1));
            end
            
            %% building P(spikes|stimuli)using kernel
            encode = nan(length(sorted),length(spksx));
            for g = 1:length(sorted)
                if ~isempty(sorted{g}) && numel(sorted{g})>=countThresh
                    tmpKernels = fitdist(sorted{g},'kernel');
                    encode(g,:) = pdf(tmpKernels,spksx);
                end
            end
            
            %% building P(stimuli|spikes) using raw data with kernel
            
            decoder = zeros(length(binrange),max(spikesAtTouch)+1);
            for k = 0:max(spikesAtTouch)
                selectspikestrials = find(spikesAtTouch == k);
                selectedtouches = varAtTouch(selectspikestrials);
                kertvar_vals = fitdist(selectedtouches,'kernel');
                decoder(:,k+1) = pdf(kertvar_vals,binrange);
                
                % varmean = mean(selectedtouches) ;
                % varstd = std(selectedtouches);
                % decoder(:,k+1) = normpdf(x,varmean,varstd);
            end
            
            %% Building P(stimuli|spikes) using bayes
            % P(stim|spks) = P(spks|stim) .* (P(stim)./P(spks))
            bloop2=parula(max(spikesAtTouch)+1);
            pstimgspks{rec} = nan(length(binrange),length(spksx));
            
            for y = 1:length(binrange)
                pstimgspks{rec}(y,:) = encode(y,:)  .*   (repmat(pdfstim(y),1,length(pdfspks)) ./  pdfspks);
            end
            
            %% Shannon Entropy
            % stimulus entropy H(S) = - sum of p(s) .* log2(p(s))
            pstimEntropy = sum(pdfstim(mintouches).*log2(pdfstim(mintouches)))*-1;
            
            %equivocation H(S|R) = -sum of p(r) .* p(stim|r) .* log2(p(stim|r));
            equivocation =   sum(repmat(pdfspks,length(mintouches),1) .*      pstimgspks{rec}(mintouches,:)      .*       log2(pstimgspks{rec}(mintouches,:))) * -1;
            
            %mutual information = reduction of uncertainty about stimulus
            %by knowing neural response. difference between H(S) and H(S|R)
            mutualInformation = pstimEntropy - equivocation;
            
            %perfect knowledge of stimuli from neuronal activity means
            %mutualInformation = pstimEntropy. If stimuli and responses are
            %independent of one another, mutual information = 0 
            
            
            
             info =  sum((repmat(pdfspks,length(mintouches),1).*pstimgspks{rec}(mintouches,:))     .*            log2(pstimgspks{rec}(mintouches,:) ./ repmat(pdfstim(mintouches)',1,length(spksx))));
            
            
            
            
            Jprob = repmat(pdfstim',1,length(spksx)) .* repmat(pdfspks,length(pdfstim),1);
            shannoninfo = Jprob(:,1).*log2(Jprob(:,1)./(pdfstim'.*pdfspks(:,1))) ;
            
            
   %%TEST model accuracy 
        [~,r] = min(abs(binrange-loo(2)));
        bayesclass(p,:) = [loo pstimgspks{rec}(:,loo(1)+1)'];
        end
           [sortedbayes, sortedBy ,binBounds]=binslin(bayesclass(:,2),bayesclass(:,3:end),'equalE',numel(V.bounds{d}),V.bounds{d}(1),V.bounds{d}(end));
            confmat=cell2mat(cellfun(@(x) mean(x,1),sortedbayes,'uniformoutput',0));
            
            %from each bin, find mean and 95% ci of the predicted theta
            %value
            ci = nan(length(sortedbayes),2);
            mid = nan(length(sortedbayes),1);
            for m = 1:length(sortedbayes) 
            [~,idx] = max(sortedbayes{m},[],2); 
            sem = std(idx)/sqrt(length(idx));
            tn = tinv([0.025 .975],length(idx)-1);
            mid(m) = mean(idx);
            ci(m,:) = mean(idx)+tn*sem;
            end

            cibounds=(ci.*(binspace))+(binrange(1)-binspace);
            midpt=(mid.*(binspace))+(binrange(1)-binspace);
            resol(rec)=nanmean(abs([binrange' - midpt]));
            resolstd(rec)=nanstd(abs([binrange' - midpt]));
            
            
            figure(580);clf;boundedline(binrange,midpt,cibounds(:,2)-midpt,'r')
            hold on;plot(binrange,binrange,'k-.')
            vals = binrange(mintouches);
            set(gca,'xlim',[vals(1) vals(end)],'ylim',[vals(1) vals(end)])
            axis square
            xlabel('Theta Actual');ylabel('Theta Predicted');
            
            
            
            
            
            figure(55);clf;imagesc(confmat')
            ylabel('Theta Predicted');xlabel('Theta Actual') 
            set(gca,'xlim',[mintouches(1) mintouches(end)],'ylim',[mintouches(1) mintouches(end)],'xtick',1:2:length(binrange),'ytick',1:2:length(binrange),'xticklabel',binrange(1:2:length(binrange)),'yticklabel',binrange(1:2:length(binrange)),'ydir','normal');
            axis square
%             %% Plotting ALL
%             %Plotting P(stimuli)
%             figure(30);subplot(5,4,rec)
%             hold on;bar(binrange,counttmp(:,1)./sum(counttmp(:,1)))
%             hold on;plot(binrange,pdfstim);
%             
%             %Plotting P(spikes)
%             figure(40);subplot(5,4,rec)
%             hold on;histogram(spikesAtTouch,'normalization','probability')
%             plot(spksx,pdfspks);
%             
%             %Plotting P(spikes|stimuli)
%             figure(50);subplot(5,4,rec)
%             bloop = parula(size(binrange,2));
%             for k = 1:length(mintouches)
%                 hold on; plot(spksx,encode(mintouches(k),:),'color',bloop(mintouches(k),:))
%             end
%             set(gca,'xtick',spksx)
%             xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
%             ylabel('Probability');
%             
%             %Plotting P(stimuli|spikes)
%             figure(60);subplot(5,4,rec)
%             bloop2=parula(max(spikesAtTouch)+1);
%             for k=1:size(decoder,2)
%                 hold on;plot(binrange,decoder(:,k),'color',bloop2(k,:))
%             end
%             set(gca,'xlim',[binrange(1) binrange(end)])
%             
%             %Plotting P(stimuli|spikes) calculated using Bayes
%             figure(90);subplot(5,4,rec)
%             hold on;
%             for p = 1:length(spksx)
%                 hold on; plot(binrange,pstimgspks{rec}(:,p),'color',bloop2(p,:))
%             end
%             set(gca,'ylim',[0 max(max(pstimgspks{rec}))])
    end
        

%     %% Aside for plotting just one cell
%     if rec == 1
%         figure(80);clf
%         subplot(2,2,1)
%         hold on;histogram(spikesAtTouch,'normalization','probability')
%         plot(spksx,pdfspks);
%         %         xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window']);
%         ylabel('P(spikes)')
%         set(gca,'xlim',[-.5 max(spikesAtTouch)+.5])
%         
%         subplot(2,2,2)
%         hold on;bar(binrange,counttmp(:,1)./sum(counttmp(:,1)))
%         plot(binrange,pdfstim);
%         %         xlabel('Whisker Angle');
%         ylabel('P(stimuli)')
%         set(gca,'xlim',[binrange(1) binrange(end)+10])
%         
%         subplot(2,2,3)
%         for k = 1:length(mintouches)
%             hold on; plot(spksx,encode(mintouches(k),:),'color',bloop(mintouches(k),:))
%         end
%         set(gca,'xtick',spksx)
%         xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
%         ylabel('P(spikes|stimuli)')
%         
%         subplot(2,2,4)
%         for p = 1:length(spksx)
%             hold on; plot(binrange,pstimgspks{rec}(:,p),'color',bloop2(p,:))
%         end
%         set(gca,'ylim',[0 max(max(pstimgspks{rec}))])
%         xlabel('Whisker Angle');ylabel('P(stimuli|spikes)')
%         
%         
%         
%     end
end

