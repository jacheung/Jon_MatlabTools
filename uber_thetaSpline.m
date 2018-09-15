%INPUT: Load uber array (L3,L5b, or L4)
Uraw=[NL5b BVL5b];

clearvars -except BVL5b NL5b U
touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));
%%
U = Uraw(selectedCells);
% U = Uraw;


modidx=cell(1,6);
[modidx{:}]=deal(zeros(1,length(U)));
valuesspline = cell(1,length(U));
valuesraw = cell(1,length(U));
window= [5 35]; %spike sum from 25ms after touch idx 


for rec = 1:length(U) 
touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
% spikesAligned = zeros(numel(touchIdx),75);
spikes = squeeze(U{rec}.R_ntk);
    for var = [1]
        theta = squeeze(U{rec}.S_ctk(var,:,:));
        thetaAtTouch = theta(touchIdx);
        thetaAtTouch = thetaAtTouch(~isnan(thetaAtTouch));
        
        %this is used to trim out touches that are way too far out
        ttranges = [mean(thetaAtTouch)-2*std(thetaAtTouch) mean(thetaAtTouch)+2*std(thetaAtTouch)];
        keepT = find([thetaAtTouch>ttranges(1) & thetaAtTouch<ttranges(2)]);
        propKept(rec) = numel(keepT)./numel(thetaAtTouch);
        
        thetaAtTouch = thetaAtTouch(keepT); 
        touchIdx = touchIdx(~isnan(thetaAtTouch)); 

        spikesAtTouch = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(thetaAtTouch),1)),2);
        
        [fo,g o] = fit(thetaAtTouch,spikesAtTouch,'smoothingspline','smoothingparam',.1);
%         figure(50);hold on; subplot(5,4,rec);
%         plot(fo,thetaAtTouch,spikesAtTouch)
        fowind=fo(linspace(min(thetaAtTouch),max(thetaAtTouch),numel(thetaAtTouch)));
        [~,indmax]=max(abs(fowind));[~,indmin]=min(abs(fowind));
        modidx{var}(rec) = (fowind(indmax)-fowind(indmin))/((fowind(indmax)+fowind(indmin)));
        
%         text(max(thetaAtTouch)*.8,max(spikesAtTouch)*.8,num2str(modidx{var}(rec)),'fontsize',20)
%         legend('off')
%         xlabel([]);ylabel([])
%         set(gca,'xlim',[min(thetaAtTouch)-5 max(thetaAtTouch)+5])
       
        if var == 1
        valuesspline{rec}(1:length(fowind),1)=fowind;
        valuesspline{rec}(1:length(fowind),2)=normalize_var(sort(thetaAtTouch),-1,1); %NORMALIZING SO ALL HAV EHTE SAME RANGE
        
        valuesraw{rec}(1:length(thetaAtTouch),2) = sort(normalize_var(thetaAtTouch,-1,1));
        [~,sortidx] = sort(normalize_var(thetaAtTouch,-1,1));
        valuesraw{rec}(1:length(thetaAtTouch),1) = spikesAtTouch(sortidx);
        end
    end
end

% a = axes;
% t1 = title(['Spike count in ' num2str(window(1)) ' to ' num2str(window(2)) ' ms post touch']);
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');

% touchNames = {'theta','phase','amp','setpoint'};
% tmp = vertcat(modidx{1},modidx{5},modidx{3},modidx{4});
% figure(5);imagesc(sortrows(tmp')');colormap(parula);hCBar=colorbar;
% set(gca,'ytick',[1:4],'yticklabel',touchNames)
% title('Modulation Index Sorted by Theta');xlabel('Sorted Cells')

%% Plotting Smoothing Splines

valsTouse = valuesraw; 
clear plotVals

%Plot raw traces
figure(11);clf
for j=1:length(valsTouse)
    hold on; plot(valsTouse{j}(:,2),valsTouse{j}(:,1));
end

%Sort splines by peak responses 
maxT = nan(length(valsTouse),1);
for k = 1:length(valsTouse)
    spks = valsTouse{k}(:,1);
    [~,pidx] = max(spks);
    maxT(k) = valsTouse{k}(pidx,2);
end
[~,sortIdx] = sort(maxT);
plotVals = valsTouse(sortIdx);

%Normalize splines of spikes 
for rec=1:length(valsTouse)
    nspks = normalize_var(plotVals{rec}(:,1),0,1);
    nspks(isnan(nspks))=0;
    plotVals{rec}(:,1) = nspks;
end

%make a heatmap
clear heat
bins = 15; 
heat=zeros(bins-1,length(valsTouse));
for rec = 1:length(valsTouse)
    [s,sb] = binslin(plotVals{rec}(:,2),plotVals{rec}(:,1),'equalE',bins,-1,1);
    heat(:,rec) = normalize_var(cellfun(@nanmean,s),0,1);
end
[~,idx ] = max(heat);
[~,heatidx] = sort(idx);


%Plot sorted and normalized smoothing spline
figure(10);clf;hold on;subplot(1,2,1)
hold on;plotVals = plotVals(heatidx);
for j=1:length(plotVals)
    plot(plotVals{j}(:,2),plotVals{j}(:,1)+j,'k');
end
set(gca,'xlim',[-1 1],'xtick',[],'ytick',[])

hold on;subplot(1,2,2)
imagesc(flipud(heat(:,heatidx)'));
set(gca,'ytick',[],'xtick',[])

% imagesc(imgaussfilt(flipud(heat(:,heatidx)'),[1 1],'padding','replicate'))


%% 
weight=cell(size(plotVals));
deg=.1;
for i=1:length(plotVals)
    for j = 1:length(plotVals{i})
    start=plotVals{i}(j,2);
    above=find(plotVals{i}(:,2)<=start+deg & plotVals{i}(:,2)>start);
    below=find(plotVals{i}(:,2)>=start-deg & plotVals{i}(:,2)<=start);
    weight{i}(j,1)=numel(vertcat(above,below));
    end
end

figure(10);clf;
 x = plotVals{1}(:,2)';
y = plotVals{1}(:,1)';
z = zeros(size(x));
col = (weight{1}/norm(weight{1},inf)*6)';  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    for k = 2:length(plotVals)
         x = plotVals{k}(:,2)';
        y = plotVals{k}(:,1)'+k;
        z = zeros(size(x));
        col = (weight{k}/norm(weight{k},inf)*6)';  % This is the color, vary with x in this case.
        hold on; surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',2);
    end
    
% for k=1:length(weight)
%     text(-8,k,num2str(max(weight{k})),'FontSize',8,'Color','black')
% end
set(gca,'ytick',[])

% colormap(jet)

% 
bonemap = bone;
bonemap = bonemap(end:-1:1,:);
colormap(bonemap);
% set(gca,'xlim',[-20 60],'xtick',-20:20:60,'visible','off')
% set(gca,'xlim',[-1 1],'xtick',[],'visible','off')
box off
