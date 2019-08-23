%col 1 = motorPosition
%col 2 = #touches at motorPosition
%col 3 = number of spikes within 30ms window of each touch at
%motorPosition
%col 4 = number of spikes from each touch 
%col 5 = theta position at touch
touchOnIdx = [find(U{rec}.S_ctk(9,:,:)==1); find(U{rec}.S_ctk(12,:,:)==1)];
touchOnIdx = sort(touchOnIdx,1);
spikesAligned = zeros(numel(touchOnIdx),31);
spikes = squeeze(U{rec}.R_ntk);

    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchOnIdx(i)+[0:30]);
    end
    
touchpos=zeros(length(U{rec}.meta.motorPosition),4);
touchpos(:,1)=U{rec}.meta.motorPosition;
touchnum=floor(touchOnIdx/4000)+1;
count=1;

for i=1:U{rec}.k
    touchpos(i,2)=sum(touchnum==i); 
    for j=count:(count-1)+max(touchpos(i,2))
        touchpos(i,3)=touchpos(i,3)+sum(spikesAligned(j,:));
    end
    count=sum(touchpos(:,2))+1;
    touchpos(i,4)=touchpos(i,3)/touchpos(i,2);
end

format shortG
touchpos=sortrows(touchpos,1);
touchpos(isnan(touchpos))=0; %set NaN in spks/touch to be zero

[sorted sortedBy binBounds]=binslin(touchpos(:,1),touchpos(:,4),'equalE',19,0,180000);

touchspks=zeros(length(sorted),1);
for i=1:length(sorted)
    touchspks(i)=sum(sorted{i})/numel(sorted{i});
end
motorpos=zeros(length(sortedBy),1);
for j=1:length(sortedBy)
    motorpos(j)=numel(sortedBy{j})/U{rec}.k;
end
figure;subplot(2,1,1)
line([5000:10000:175000],touchspks);
hold on
scatter([5000:10000:175000],touchspks);
set(gca,'xtick',[0:10000:180000]);
xlabel('Motor Position (um)')
ylabel('spks/touch')
subplot(2,1,2)
bar([5000:10000:175000],motorpos)
set(gca,'xtick',[0:10000:180000]);
xlabel('Motor Position (um)')
ylabel('Proportion of Trials')

%% Theta at touch

theta=U{rec}.S_ctk(1,:,:);
touchtheta(:,1)=touchOnIdx;
touchtheta(:,2)=theta(touchOnIdx);
spikesAligned = zeros(numel(touchOnIdx),31);
window = [5 10 20 30];

    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchOnIdx(i)+[0:window(4)]);
        touchtheta(i,3)=sum(spikesAligned(i,:));
    end
    
sortrows(touchtheta,2);
[sorted sortedBy binBounds]=binslin(touchtheta(:,2),touchtheta(:,3),'equalE',71,-20,50);

for i=1:length(sorted)
    spkmean(i)=sum(sorted{i})/numel(sorted{i});
end
spkmean(isnan(spkmean))=0;
figure;bar([-19.5:1:49.5],spkmean)
title('theta at touch vs spks')
ylabel('spks/touch at theta')
xlabel('theta at touch')



