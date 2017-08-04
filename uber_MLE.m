touchOnIdx = [find(U{rec}.S_ctk(9,:,:)==1); find(U{rec}.S_ctk(12,:,:)==1)];
touchOffIdx = [find(U{rec}.S_ctk(10,:,:)==1); find(U{rec}.S_ctk(13,:,:)==1)];
spikes = squeeze(U{rec}.R_ntk);
touchOnIdx = touchOnIdx(touchOnIdx<(numel(spikes)-50));
touchOffIdx = touchOffIdx(1:length(touchOnIdx));
touchOnIdx = sort(touchOnIdx,1);
touchOffIdx = sort(touchOffIdx,1);
touchrange=horzcat(touchOnIdx,touchOffIdx);

%align spikes all around touches (-25ms:50ms after touch)
spikesAligned = zeros(numel(touchOnIdx),76);
spikes = squeeze(U{rec}.R_ntk);
for i = 1:size(spikesAligned,1)
    spikesAligned(i,:) = spikes(touchOnIdx(i)+[-25:50]);
end

%% Normal dist fit
%need to feed it time of when spikes occurred
inds=find(spikesAligned==1); %find index where all spikes occurred w/in touch window
[~, col] = ind2sub(size(spikesAligned),inds); %gets column associated with index (ie. time in ms of where spike occurred)
col(col<50)=NaN;%set all spikes before 0 to NaN and not be counted in distribution. Do this so we can get a distribution that truly matches touch response
col(col>75)=NaN;
dist=fitdist(col,'Normal'); %we just want to find a good fit for the touch response. Doesnt have to be poisson here as we use poisson to dictate spike probability in the GLM further down. 
dist.mu=dist.mu-50; %to account for index starting at 0 instead of -50

figure;
bar(-25:50,sum(spikesAligned)/numel(spikesAligned),'k')
title('First Touch Onset')
xlabel('Time from touch (ms)')
ylabel('Spikes per touch')
hold on
plot([-25:50],pdf(dist,[-25:50])/40,'LineWidth',2,'color','g')
%% Touch spike design matrix 
epoch=spikesAligned(:); 

%% Theta stim design matrix 
samp=U{rec}.S_ctk(1,touchOnIdx);%find theta at all touches 
psamp=U{rec}.S_ctk(5,touchOnIdx);
ampsamp=U{rec}.S_ctk(3,touchOnIdx);
spsamp=U{rec}.S_ctk(4,touchOnIdx);
thetaAligned=zeros(length(samp),length(-25:50));
phaseAligned=zeros(length(psamp),length(-25:50));
ampAligned=zeros(length(ampsamp),length(-25:50));
spAligned=zeros(length(spsamp),length(-25:50));


for i=1:length(samp)
     [~,I] = min(abs(samp(i)-thetatouchspks(:,1)));
     thetaAligned(i,:) = smooth(thetatouchspks(I,2:77));
     [~,I] = min(abs(psamp(i)-phasetouchspks(:,1)));
     phaseAligned(i,:) = smooth(phasetouchspks(I,2:77));
     [~,I] = min(abs(ampsamp(i)-amptouchspks(:,1)));
     ampAligned(i,:) = smooth(amptouchspks(I,2:77));
     [~,I] = min(abs(spsamp(i)-sptouchspks(:,1)));
     spAligned(i,:) = smooth(sptouchspks(I,2:77));
end

thetaAligned=reshape(thetaAligned',[],1);
phaseAligned=reshape(phaseAligned',[],1);
ampAligned=reshape(ampAligned',[],1);
spAligned=reshape(spAligned',[],1);

% spike history stimulus
% not sure how to create the stimulus for this yet. 
% GLM BUILD 
[gnew,dev,stats] = glmfit([phaseAligned ampAligned spAligned],thetaAligned,'poisson');
figure; plot(1:750,exp(gnew(1)+gnew(2)*phaseAligned(1:750)+gnew(3)*ampAligned(1:750)+gnew(4)*spAligned(1:750)))
hold on;plot(1:750,thetaAligned(1:750),'r')
%% Touch stim design matrix 
hypparam=dist.mu
tft=repmat([(-50:49)'], [length(touchOnIdx), 1]); %time from touch (50ms before to 50ms after)
newtft=exp(-(tft-hypparam).^2/(2*dist.sigma^2));%tft into gaussian distrubtion w/ mean at 10ms from touch

%% GLM BUILD
[gnew,dev,stats] = glmfit([newtft thetaAligned phaseAligned ampAligned spAligned],epoch,'poisson');%glm fit w/ modeled Gaussian dist.
%figure; plot(-50:49,exp(gnew(1)+gnew(2)*(newtft(1:100)))/0.001+gnew(3)*(thetaAligned(1:100))/0.001)
figure; plot(1:length(epoch),exp(gnew(1)+(gnew(2)*newtft)/.01+(gnew(3)*thetaAligned)/.01+(gnew(4)*phaseAligned)/.01+(gnew(5)*ampAligned)/.01+(gnew(6)*spAligned)/.01))
hold on; plot(1:length(epoch),epoch,'r')
title('GLM with Gaussian Fit')
xlabel('Time from Touch Onset (ms)')
ylabel('Firing Rate (spks/s)')
%% GLM with JUST TOUCH 
[gnew2,dev,stats] = glmfit(newtft,epoch,'poisson');%glm fit w/ modeled Gaussian dist.
figure; plot(1:53600,exp(gnew2(1)+gnew(2)*(newtft)));
hold on; plot(1:53600,epoch,'r')

%% AIC for model selection. Lower AIC = better model 
%use exp(a0)*Delta to calculate FR/ms as the results are in log. exp gets rid of
delt=.001;%for ms

ratefunc=(gnew2(1)+gnew2(2)*newtft);
params=length(gnew2);
llfunc=sum(epoch.*(ratefunc)-delt*exp(ratefunc)); %log likelihood function
aic=(-2*llfunc)+(2*params); % aic = -2loglikelihod + 2x(params)

allratefunc=(gnew(1)+gnew(2)*newtft+gnew(3)*thetaAligned+gnew(4)*phaseAligned+gnew(5)*ampAligned+gnew(6)*spAligned);
allparams=length(gnew);
allfunc=sum(epoch.*(allratefunc)-delt*exp(allratefunc));
allaic=(-2*allfunc)+(2*allparams);


%% fminunc
%func=@(x)