%INPUT: Load uber array (L3,L5b, or L4)

modidx=cell(1,6);
[modidx{:}]=deal(zeros(1,length(U)));
values = cell(1,length(U));

window= [8 50]; %spike sum from 25ms after touch idx 

figure(50);clf;
for rec = 1:length(U) 
touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
% spikesAligned = zeros(numel(touchIdx),75);
spikes = squeeze(U{rec}.R_ntk);
    for var = [1]
        theta = squeeze(U{rec}.S_ctk(var,:,:));
        thetaAtTouch = theta(touchIdx);
        thetaAtTouch = thetaAtTouch(~isnan(thetaAtTouch));
        touchIdx = touchIdx(~isnan(thetaAtTouch)); 

        spikesAtTouch = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(thetaAtTouch),1)),2);
        
        [fo,g o] = fit(thetaAtTouch,spikesAtTouch,'smoothingspline','smoothingparam',.1);
%         figure(50);hold on; subplot(5,4,rec);
%         plot(fo,thetaAtTouch,spikesAtTouch)
        fowind=fo(linspace(min(thetaAtTouch),max(thetaAtTouch),numel(thetaAtTouch)));
        [~,indmax]=max(abs(fowind));[~,indmin]=min(abs(fowind));
        modidx{var}(rec) = (fowind(indmax)-fowind(indmin))/((fowind(indmax)+fowind(indmin)));
        
        text(max(thetaAtTouch)*.8,max(spikesAtTouch)*.8,num2str(modidx{var}(rec)),'fontsize',20)
        legend('off')
        xlabel([]);ylabel([])
        set(gca,'xlim',[min(thetaAtTouch)-5 max(thetaAtTouch)+5])
       
        if var == 1
        values{rec}(1:length(fowind),1)=fowind;values{rec}(1:length(fowind),2)=sort(thetaAtTouch);
        end
    end
end

a = axes;
t1 = title(['Spike count in ' num2str(window(1)) ' to ' num2str(window(2)) ' ms post touch']);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');

% touchNames = {'theta','phase','amp','setpoint'};
% tmp = vertcat(modidx{1},modidx{5},modidx{3},modidx{4});
% figure(5);imagesc(sortrows(tmp')');colormap(parula);hCBar=colorbar;
% set(gca,'ytick',[1:4],'yticklabel',touchNames)
% title('Modulation Index Sorted by Theta');xlabel('Sorted Cells')

%% Plotting Smoothing Splines
%Plot raw traces
figure(10);plot(values{1}(:,2),values{1}(:,1));
hold on;
for j=2:length(values)
    plot(values{j}(:,2),values{j}(:,1));
end

%Sort splines by peak responses 
test2=values;
test3=cell(1,length(U));
peaktheta=zeros(2,length(U));

for rec=1:length(U)
    if numel(find(test2{rec}(:,1)==max(test2{rec}(:,1))))>1
        peaktheta(1,rec)=NaN;
    else peaktheta(1,rec) = test2{rec}(find(test2{rec}(:,1)==max(test2{rec}(:,1))),2);
    end
end
peaktheta(2,:)=1:length(peaktheta);
sortedpeak=sortrows(peaktheta',1);
for rec=1:length(U)
    test3{rec}=test2{sortedpeak(rec,2)};
end

%Normalize splines and then plot 
for rec=1:length(U)
        tmp7=test3{rec}(:,1)-min(test3{rec}(:,1));
        tmp8=tmp7/((max(test3{rec}(:,1))-min(test3{rec}(:,1))));
        if isnan(tmp8) %added this in case there were no spikes at all as 0/0 = NaN
            tmp8=0; %this ensures that even those with no spikes gets plotted
        end
        test3{rec}(:,1)=tmp8+(rec);%+ rec so each smoothed spline is plotted on a different y axis
end

%Plot sorted and normalized smoothing spline
figure(9);clf;plot(test3{1}(:,2),test3{1}(:,1));
hold on;
for j=2:length(test3)
    plot(test3{j}(:,2),test3{j}(:,1));
end
set(gca,'xlim',[-30 60])
%%
weight=cell(size(test3));
deg=1;
for i=1:length(test3)
    for j = 1:length(test3{i})
    start=test3{i}(j,2);
    above=find(test3{i}(:,2)<=start+deg & test3{i}(:,2)>start);
    below=find(test3{i}(:,2)>=start-deg & test3{i}(:,2)<=start);
    weight{i}(j,1)=numel(vertcat(above,below));
    end
end

figure(10);clf;
x = test3{1}(:,2)';
y = test3{1}(:,1)';
z = zeros(size(x));
col = (weight{1}/norm(weight{1},inf)*6)';  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    for k = 2:length(test3)
        x = test3{k}(:,2)';
        y = test3{k}(:,1)';
        z = zeros(size(x));
        col = (weight{k}/norm(weight{k},inf)*6)';  % This is the color, vary with x in this case.
        hold on; surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',2);
    end
    
% for k=1:length(weight)
%     text(-8,k,num2str(max(weight{k})),'FontSize',8,'Color','black')
% end

set(gca,'ytick',[])
bonemap = bone;
bonemap = bonemap(end:-1:1,:);
colormap(bonemap);
set(gca,'xlim',[-20 60],'xtick',-20:20:60,'visible','off')
box off
