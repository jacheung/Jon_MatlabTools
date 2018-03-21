
%% Build for a ghetto decoder.
%P(spks|stimuli) given as just a poisson distribution with lambda of the
%mean of all spks generated within the window for that angle bin
%P(stimuli|spks) given as just a normal distribution with mean of all
%stimuli that gave x number of spks.
%Decoder samples spikes for each bin based on P(spks|stimuli) and finds the
%P(bin|spks) using the normal distributions.




window = [8 20];
countThresh = 10;
trials2sample = 100;
neuronsample = 200;

figure(30);clf
figure(40);clf;
figure(50);clf;
figure(60);clf;
figure(80);clf;
figure(100);clf;
clear decodeAcc
clear pstimgspks


for rec = [1 3 4 5 9 10 11 12 14]
    
    
    
   
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
    V.bounds = [{[-40:2.5:80]} {[-9750:125:9750]} {[-99.5:5:99.5]} {[-99.5:5:99.5]} {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
    
    for d = 1
        var = squeeze(U{rec}.S_ctk(d,:,:));
        varAtTouch = var(touchIdx);
        varAtTouch = varAtTouch(~isnan(varAtTouch));
        touchIdx = touchIdx(~isnan(varAtTouch));
        spikesAtTouch = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(varAtTouch),1)),2);
        
        [sorted, sortedBy ,binBounds]=binslin(varAtTouch,spikesAtTouch,'equalE',numel(V.bounds{d})+1,V.bounds{d}(1),V.bounds{d}(end));
        binrange = V.bounds{d};
        trims=[binrange' cell2mat(cellfun(@mean,sorted,'Uniformoutput',0)) cell2mat(cellfun(@std,sorted,'Uniformoutput',0))];
        counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
        mintouches=find(counttmp(:,1)>countThresh);
        
        %% norm distribution for building P(spks)
        spksx = [0:1:ceil(max(spikesAtTouch))];
%         pdfspks = normpdf(spksx,mean(spikesAtTouch),std(spikesAtTouch));
        
                ker_vals = fitdist(spikesAtTouch,'kernel');
                pdfspks = pdf(ker_vals,spksx);
        
        
        figure(40);subplot(5,4,rec)
        hold on;histogram(spikesAtTouch,'normalization','probability')
        plot(spksx,pdfspks);
        
        
        
        %% norm distribution for building P(stimuli)
        stims = [];
        for y = 1:length(binrange)
            stims = [stims ; repmat(binrange(y),counttmp(y,1),1)];
        end
        %kernel distribution 
          kervar_vals = fitdist(varAtTouch,'kernel');
          pdfstim = pdf(kervar_vals,binrange);
          
          %normal distribution
%         pdfstim = normpdf(binrange,mean(varAtTouch),std(varAtTouch));
        figure(30);subplot(5,4,rec)
        hold on;bar(binrange,counttmp(:,1)./sum(counttmp(:,1)))
        hold on;plot(binrange,pdfstim);
        
        
        
        %% poisson distribution for building P(spikes|stimuli)
%         x = [0:1:ceil(max(spikesAtTouch))];
%         lambda = trims(:,2);
%         encode = nan(length(lambda),length(x));
%         for g = 1:length(lambda)
%             encode(g,:) = poisspdf(x,lambda(g));
%         end
%         sampledSpikes = poissrnd(repmat(lambda,1,trials2sample));
%         
%         figure(50);subplot(5,4,rec)
%         bloop = parula(size(binrange,2));
%         for k = 1:length(mintouches)
%             hold on; plot(x,encode(mintouches(k),:),'color',bloop(mintouches(k),:))
%         end
%         set(gca,'xtick',x)
%         xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
%         ylabel('Probability')
        
        %OR 
        % kernel distribution for building P(spikes|stimuli)     
        spksx = [0:1:ceil(max(spikesAtTouch))];
        
        encode = nan(length(sorted),length(spksx));
        for g = 1:length(sorted)
            if ~isempty(sorted{g}) && numel(sorted{g})>=countThresh
            tmpKernels = fitdist(sorted{g},'kernel');
            encode(g,:) = pdf(tmpKernels,spksx);
            end
        end
        

        figure(50);subplot(5,4,rec)
        bloop = parula(size(binrange,2));
        for k = 1:length(mintouches)
            hold on; plot(spksx,encode(mintouches(k),:),'color',bloop(mintouches(k),:))
        end
        set(gca,'xtick',spksx)
        xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
        ylabel('Probability')
        
        
        %% norm distribution for building P(stimuli|spikes)
        %decoder is a m x n matrix (m = number of possible spikes 0:x and n
        %= binranges)
        x = [binrange];
%         kl=zeros(max(spikesAtTouch)+1,1);
        
        bloop2=parula(max(spikesAtTouch)+1);
        decoder = zeros(max(spikesAtTouch)+1,length(x));
        for k = 0:max(spikesAtTouch)
            selectspikestrials = find(spikesAtTouch == k);
            selectedtouches = varAtTouch(selectspikestrials);
            
            %norm distribution 
%              varmean = mean(selectedtouches) ;
%             varstd = std(selectedtouches);
%             decoder(k+1,:) = normpdf(x,varmean,varstd);
            
            %kernel distribution 
            kertvar_vals = fitdist(selectedtouches,'kernel');
            decoder(k+1,:) = pdf(kertvar_vals,binrange);

            figure(60);subplot(5,4,rec)
            hold on;plot(x,decoder(k+1,:),'color',bloop2(k+1,:))
%             kl(k+1) = filex_kldiv(binrange,pdfstim,decoder(k+1,:));
        end
        
        %plotting kuelback lieblier distribution for rate
%         figure(100)
%         hold on; plot(0:1/(numel(kl)-1):1,kl)
        
        
%         decodedProb = nan(size(sampledSpikes));
%         for u = 1:size(sampledSpikes,1)
%             currSpikes = sampledSpikes(u,:);
%             if ~isnan(currSpikes)
%                 for g = 1:numel(currSpikes)
%                     if currSpikes(g)+1<size(decoder,1)
%                         decodedProb(u,g) = decoder(currSpikes(g)+1,u);
%                     end
%                 end
%             end
%         end
    end
    
%     decodeAcc(:,rec) = nanmean(decodedProb,2);
    
    %% BAYES OPTIMAL CLASSIFIER 
    % P(stim|spks) = P(spks|stim) .* (P(stim)./P(spks))
    pstimgspks{rec} = nan(length(binrange),length(spksx));
    
    for y = 1:length(binrange)
    pstimgspks{rec}(y,:) = encode(y,:)  .*   (repmat(pdfstim(y),1,length(pdfspks)) ./  pdfspks); %P(spikes|stimuli)
    end
    
    figure(110);subplot(5,4,rec)
    hold on; 
    for p = 1:length(spksx)
        hold on; plot(1:length(binrange),smooth(pstimgspks{rec}(:,p)),'color',bloop2(p,:))
    end
    set(gca,'ylim',[0 max(max(pstimgspks{rec}))])
    %% Aside for plotting just one cell
    if rec == 1
        figure(80);clf
        subplot(2,2,1)
        hold on;histogram(spikesAtTouch,'normalization','probability')
        plot(spksx,pdfspks);
        xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window']); ylabel('P(spikes)')
    
        
        subplot(2,2,2)
        hold on;bar(binrange,counttmp(:,1)./sum(counttmp(:,1)))
        plot(binrange,pdfstim);
        xlabel('Whisker Angle');ylabel('P(stimuli)')
        set(gca,'xlim',[binrange(1) binrange(end)+10])
        
        subplot(2,2,3)
        for k = 1:length(mintouches)
            hold on; plot(spksx,encode(mintouches(k),:),'color',bloop(mintouches(k),:))
        end
        set(gca,'xtick',spksx)
        xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
        ylabel('P(spikes|stimuli)')
        
        subplot(2,2,4)
%         for f = 0:max(spikesAtTouch)
%             hold on;plot(binrange,decoder(f+1,:),'color',bloop2(f+1,:))
%         end
        xlabel('Whisker Angle');ylabel('P(stimuli|spikes)') 
        
        for p = 1:length(binrange)
        hold on; plot(spksx,pstimgspks{rec}(p,:),'color',bloop(p,:))
        end
        
        
        figure(90);clf
        subplot(2,1,1)
        plot(binrange,pdfstim,'k','linewidth',5);

        for f = 0:max(spikesAtTouch)
            hold on;plot(binrange,decoder(f+1,:),'color',bloop2(f+1,:))
        end
        legend({'P(stim)' num2str(1:numel(kl))})
        subplot(2,1,2)
            plot(0:1/(numel(kl)-1):1,kl)
            xlabel('normalized rate') 
            ylabel('KL value')      
        
    end
end


%% POPULATION DECODER... NOT VERY GOOD. QUITE GHETTO AND ACTUALLY INCORRECT 180221
popdecoder = zeros(size(decodeAcc,1),neuronsample);
for g = 1:size(decodeAcc,1)
    probs = decodeAcc(g,~isnan(decodeAcc(g,:)));
    
    if ~isempty(probs)
        tmp=datasample(probs,neuronsample);
        serialsums = zeros(length(tmp),length(tmp)*2);
        for k = 0:length(tmp)
            serialsums(k+1,(k+1:k+length(tmp))) = tmp;
        end
        
        popdecoder(g,:)= nansum(serialsums(:,1:length(tmp)));
    end
end
figure(680);clf
for t = 1:size(popdecoder,1)
    hold on; plot(1:size(popdecoder,2),popdecoder(t,:),'color',bloop(t,:));
end
legend(num2str(binrange'))
set(gca,'ylim',[0 1])
