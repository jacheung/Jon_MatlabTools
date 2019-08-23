
%% find all indices where touch occurs and create window around it. 
touchOnIdx = [find(U{rec}.S_ctk(9,:,:)==1); find(U{rec}.S_ctk(12,:,:)==1)];
touchOffIdx = [find(U{rec}.S_ctk(10,:,:)==1); find(U{rec}.S_ctk(13,:,:)==1)];
spikes = squeeze(U{rec}.R_ntk);
touchOnIdx = touchOnIdx(touchOnIdx<(numel(spikes)-50)); %limits touches to those within recording window
touchOffIdx = touchOffIdx(1:length(touchOnIdx));
touchOnIdx = sort(touchOnIdx,1);
touchOffIdx = sort(touchOffIdx,1);
touchrange=horzcat(touchOnIdx,touchOffIdx); %all touch windows

%align spikes all around touches (-25ms:50ms after touch)
spikesAligned = zeros(numel(touchOnIdx),76);
spikes = squeeze(U{rec}.R_ntk);
for i = 1:size(spikesAligned,1)
    spikesAligned(i,:) = spikes(touchOnIdx(i)+[-25:50]);
end

%% Identify Dependencies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atTouchvar = zeros(length(touchOnIdx),3);

for i = 1:length(touchOnIdx)
    atTouchVar(i,1) = U{rec}.S_ctk(3,touchOnIdx(i)); %amplitude
    atTouchVar(i,2) = U{rec}.S_ctk(4,touchOnIdx(i)); %setpoint
    atTouchVar(i,3) = U{rec}.S_ctk(5,touchOnIdx(i)); %phase
    atTouchVar(i,4) = U{rec}.S_ctk(1,touchOnIdx(i)); %theta
end

%relevancy check (running each IV:DV)
[rho,pval]=corr(atTouchVar,'type','pearson');

figure(10);clf;subplot(2,3,1);scatter(atTouchVar(:,1),atTouchVar(:,4)) %amp x theta
    ylabel('theta');xlabel('amplitude');lsline
    text(20,6,num2str(rho(1,4)),'FontSize',8,'Color','black')
    text(20,4.5,['p=' num2str(pval(1,4))],'FontSize',8,'Color','black')
subplot(2,3,2);scatter(atTouchVar(:,2),atTouchVar(:,4))%setpoint x theta
    title('Relevancy Check');xlabel('setpoint');lsline
    text(20,6,num2str(rho(2,4)),'FontSize',8,'Color','black')
    text(20,4.5,['p=' num2str(pval(2,4))],'FontSize',8,'Color','black')
subplot(2,3,3);scatter(atTouchVar(:,3),atTouchVar(:,4)) %phase x theta
    xlabel('phase');lsline
    text(2,5,num2str(rho(3,4)),'FontSize',8,'Color','black')
    text(2,3.5,['p=' num2str(pval(3,4))],'FontSize',8,'Color','black')
% poor dependency of calculating phase but this may be due to the fact that
% the relationship is NOT linear

%multicolinearity check (running each IV:IV)
subplot(2,3,4);scatter(atTouchVar(:,1),atTouchVar(:,2)) %amp x setpoint
    xlabel('amplitude');ylabel('setpoint');lsline
    text(20,6,num2str(rho(1,2)),'FontSize',8,'Color','black')
    text(20,4.5,['p=' num2str(pval(1,2))],'FontSize',8,'Color','black')
subplot(2,3,5);scatter(atTouchVar(:,2),atTouchVar(:,3))%setpoint x phase
    xlabel('setpoint'),ylabel('phase');title('Multicollinearity');lsline
    text(40,-3,num2str(rho(2,3)),'FontSize',8,'Color','black')
    text(40,-3.25,['p=' num2str(pval(2,3))],'FontSize',8,'Color','black')
subplot(2,3,6);scatter(atTouchVar(:,1),atTouchVar(:,3)) %amp x phase
    xlabel('amplitude'),ylabel('phase');lsline
    text(20,-3,num2str(rho(1,3)),'FontSize',8,'Color','black')
    text(20,-3.25,['p=' num2str(pval(1,3))],'FontSize',8,'Color','black')

%figure(6);scatter3(atTouchVar(:,1),atTouchVar(:,2),atTouchVar(:,3))


%% GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualize heatmap in 3D
prms=repmat(thetatouchspks(:,1),1,76);
prms2=repmat(-25:50,11,1);
z=thetatouchspks(:,2:77);
figure(11);surf(prms,prms2,z); %3D plot for heatmap

xlabel('Theta at Touch'); ylabel('Time from touch (ms)'); zlabel('Change in FR') 
title('Theta Heatmap in 3D')
%% stimulus matrix 
samp=U{rec}.S_ctk(1,touchOnIdx);%find theta at all touches 
psamp=U{rec}.S_ctk(5,touchOnIdx);
ampsamp=U{rec}.S_ctk(3,touchOnIdx);
spsamp=U{rec}.S_ctk(4,touchOnIdx);
thetaAligned=zeros(length(samp),length(-25:50));
phaseAligned=zeros(length(psamp),length(-25:50));
ampAligned=zeros(length(ampsamp),length(-25:50));
spAligned=zeros(length(spsamp),length(-25:50));


for i=1:length(samp) %find bin in heatmap closest to value at touch 
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

% spike history stimulus: not sure how to create the stimulus for this yet. 
% GLM BUILD 
[gnew,dev,stats] = glmfit([phaseAligned ampAligned spAligned],thetaAligned,'poisson');
figure(31); plot(1:750,exp(gnew(1)+gnew(2)*phaseAligned(1:750)+gnew(3)*ampAligned(1:750)+gnew(4)*spAligned(1:750)))
hold on;plot(1:750,thetaAligned(1:750),'r')

text(500, .3,['Baseline FR: ' num2str(exp(gnew(1))*1000) 'spks/s'],'FontSize',8,'Color','black')
text(500, .285,['Phase Wt: ' num2str(gnew(2))],'FontSize',8,'Color','black')
text(500, .270,['Amp Wt: ' num2str(gnew(3))],'FontSize',8,'Color','black')
text(500, .255,['Setpoint Wt: ' num2str(gnew(4))],'FontSize',8,'Color','black')
xlabel('Time (ms); each trial = 75ms'), ylabel('Change in FR')
title('GLM output, Blue = raw data, Red = GLM predicted')

