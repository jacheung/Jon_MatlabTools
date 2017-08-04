[mask]= assist_touchmasks(U{rec});

spikes = squeeze(U{rec}.R_ntk);
firsttouchIdx = [find(U{rec}.S_ctk(9,:,:)==1)];
subsetouchIdx = [find(U{rec}.S_ctk(12,:,:)==1)];
firsttouchOffIdx = [find(U{rec}.S_ctk(10,:,:)==1)];
touchOnIdx = [find(U{rec}.S_ctk(9,:,:)==1); find(U{rec}.S_ctk(12,:,:)==1)];
touchOffIdx = [find(U{rec}.S_ctk(10,:,:)==1); find(U{rec}.S_ctk(13,:,:)==1)];
if touchOffIdx<(numel(spikes)-70) %ensures that offset is within index boundaries
    touchOffIdx = touchOffIdx+70;
end
motorpos=U{rec}.meta.motorPosition;
%%% Phase at touch vs subsequent touches

FTphase=zeros(1,length(firsttouchIdx));
for i = 1:length(firsttouchIdx)
    FTphase(i)=U{rec}.S_ctk(5,firsttouchIdx(i));
end
STphase=zeros(1,length(subsetouchIdx));
for i = 1:length(subsetouchIdx)
    STphase(i)=U{rec}.S_ctk(5,subsetouchIdx(i));
end
figure(99)
scatter(1:length(FTphase),FTphase,'b')
hold on;scatter(length(FTphase)+1:1:length(STphase)+length(FTphase),STphase,'r')


%%% Change in Pretouch velocity 5ms around touch
velDelt=zeros(length(touchOnIdx),11);
for i = 1:length(touchOnIdx)
    velDelt(i,1:11)=U{rec}.S_ctk(2,touchOnIdx(i)-5:touchOnIdx(i)+5);
end
mean(velDelt);

%sort velocity change x position of pole 
TNumatTouch=ceil(touchOnIdx/4000);
MotorPosatTouch = motorpos(TNumatTouch);
normvel=zeros(length(velDelt),12);
normvel(:,12)=MotorPosatTouch;

velmax=repmat(max(velDelt(:,1:11),[],2),1,11);
velmin=repmat(min(velDelt(:,1:11),[],2),1,11);
normvel(:,1:11)=(velDelt-velmin)./(velmax-velmin);
newnormvel=sortrows(normvel,12);

figure(22);
for i = 500:525
    plot(1:11,newnormvel(i,1:11)+i-1);hold on
end