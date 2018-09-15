
%%
clear pop
for rec=1:length(U)
    countThresh = 10; %min touch in each bin to be considered
    gaussFilt = [1]; %sigma filter for imgaussfilt (can be single value for filtering across rows and columns or two values [x y] for filter of row/column
    window = [-50:50]; %ms to look around touch
    
    [varspikes, preDvarspikes,postDvarspikes] = assist_varAtTouch(U{rec},window);
    % First 6 columns will be values for the variables
    % 1) THETA
    % 2) PRE TOUCH VELOCITY
    % 3) AMP
    % 4) SETPOINT
    % 5) PHASE
    % 6) MAX KAPPA
    % Last columns will be the spikes around your given window
    
    
    % Run this line below if you want to elim all protraction touches and
    % just look at ret touches.
    %     varspikes(varspikes(:,5)<0,:)=[];
    
    
    
    %     comp={varspikes, preDvarspikes, postDvarspikes};
    comp={varspikes};
    for g = 1:length(comp)
        %%
        %Plot theta at touch
        fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
        V.bounds = [{[-100:5:100]} {[-9750:125:9750]} {[-99.5:5:99.5]} {[-99.5:5:99.5]} {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
        for d = [1 2 3 4 5 6] %for variables 1:6
            
            if d == 2
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-10000,10000);
            elseif d == 5
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalX',numel(V.bounds{d})+1);
            elseif d == 6
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-1,1);
            else
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-100,100);
            end
            
            binrange = V.bounds{d};
            
            % Trimming unused ranges
            trims=[binrange' cell2mat(cellfun(@mean,sortedBy,'Uniformoutput',0))];
            indexHasTheta = ~isnan(trims(:,2));
            trims = trims(indexHasTheta, :);
            counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
            
            % Populating fields with values
            V.(fields{d}).counts =counttmp(indexHasTheta==1,1);
            V.(fields{d}).range = trims(:,1);
            V.(fields{d}).spikes=cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
            V.(fields{d}).spikes=V.(fields{d}).spikes(indexHasTheta==1,:);
            V.(fields{d}).raw = sorted(indexHasTheta==1);
            
            %Trimming bins with touchcounts below thresholds
            mintouches=find(V.(fields{d}).counts<countThresh);
            V.(fields{d}).counts(mintouches,:)=[];
            V.(fields{d}).spikes(mintouches,:)=[];
            V.(fields{d}).range(mintouches,:)=[];
            V.(fields{d}).raw(mintouches,:) = [];
            
            %Plotting features
            %             figure(30+rec);subplot(2,3,d);
            %
            %
            %
            % %             imagesc(imgaussfilt(V.(fields{d}).spikes,gaussFilt,'padding','replicate')); %GAUSS PLOTTING
            %             imagesc(V.(fields{d}).spikes) %OLD PLOTTING
            %
            %             colormap(gca,parula);
            %             set(gca,'Ydir','normal','ytick',(1:length(V.(fields{d}).range)),'yticklabel',[V.(fields{d}).range],...
            %                 'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
            %             for k=1:size(V.(fields{d}).counts,1)
            %                 text(20,k,num2str(V.(fields{d}).counts(k)),'FontSize',8,'Color','white')
            %             end
            %             hold on;plot([sum(window<0) sum(window<0)],[length(V.(fields{d}).range) 0],'w:')
            %             axis('square')
            %             xlabel('time from touch onset (ms)')
            %             title([fields{d}])
            %
            
            pop{rec}=V;
            
            
        end
        
    end
end

%% POPULATION PLOTTING

% PLOT TOUCH CELLS+

fields = [    {'theta'}  ];

% touchCells = touchCell(U,1.5);
touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);

[rc] = numSubplots(length(selectedCells));
figure(333);clf
figure(322);clf
clear popzscore

plotrow = rc(1);
plotcolumn = rc(2);
for d = 1:length(selectedCells)
    
    currCell = selectedCells(d);
    figure(333);subplot(plotrow,plotcolumn,d);
    imagesc(pop{currCell}.(fields{1}).spikes)
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{currCell}.(fields{1}).range)),'yticklabel',[pop{currCell}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{currCell}.(fields{1}).counts,1)
        text(20,k,num2str(pop{currCell}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{currCell}.(fields{1}).range) 0],'w:')
    
    
    figure(322);subplot(plotrow,plotcolumn,d);
    imagesc(imgaussfilt(pop{currCell}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{currCell}.(fields{1}).range)),'yticklabel',[pop{currCell}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{currCell}.(fields{1}).counts,1)
        text(20,k,num2str(pop{currCell}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{currCell}.(fields{1}).range) 0],'w:')
    
    %ZSCORREEEEE
    thetabounds = [-100:5:100]; 
    counts = pop{currCell}.theta.counts;
    ranges = pop{currCell}.theta.range ;
    thetavect = [];
    for k = 1:length(ranges)
        thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
    end
    
    allspikes = cell2mat(pop{currCell}.theta.raw);
    meanbase = mean(mean(allspikes(:,1:51),2));
    stdbase = std(mean(allspikes(:,1:51),2));
    postresponses = allspikes(:,55:80);
    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    uniquethetas = unique(thetavect);
    
    [sorted,sortedBy,~] = binslin(thetavect,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
    
    zraw = nan(length(sorted),1000);
    for b = 1:length(sorted)
        currvals = sorted{b};
        if ~isempty(sorted{b})
            zraw(b,1:length(currvals)) = currvals';
        end
    end
    

    [p,~,stats]=anova1(zraw',[],'off');
    vals =multcompare(stats);
    
    popzscore{d} = cellfun(@mean,sorted);
    otune(d)=p;
    
    
   tmp =  vals(vals(:,end)<.01,:);
    figure;scatter(tmp(:,1),tmp(:,2))
    set(gca,'xlim',[0 15],'ylim',[0 15])
end
%Plotting Zscore

ocells = find(otune<.01);

zs = cell2mat(popzscore(ocells))';
xticks = thetabounds(:,~nansum(zs)==0);
zs(:,nansum(zs)==0)=[];
[~,idx] = max(zs,[],2);

[~,peakidx] = sort(idx);


 figure(480);clf;imagesc(zs(peakidx,:)) %sorted by peak idx
% figure(480);clf;imagesc(zs) %unsorted
h=colorbar;
colormap(redbluecmap)
caxis([-3 3])
set(gca,'xtick',[1:2:size(zs,2)],'xticklabel',xticks(1:2:end),'ytick',[])
xlabel('Whisker angle at touch') 
    
    

%% PLOT ALL
[rc] = numSubplots(length(pop));

plotrow = rc(1);
plotcolumn = rc(2);

for d = 1:length(pop)
    
    
    figure(30);subplot(plotrow,plotcolumn,d);
    %       imagesc(imgaussfilt(pop{d}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    imagesc(pop{d}.(fields{1}).spikes)
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{d}.(fields{1}).range)),'yticklabel',[pop{d}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{d}.(fields{1}).counts,1)
        text(20,k,num2str(pop{d}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{d}.(fields{1}).range) 0],'w:')
    
    
    figure(31);subplot(plotrow,plotcolumn,d);
    imagesc(imgaussfilt(pop{d}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    %      imagesc(pop{d}.(fields{1}).spikes)
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{d}.(fields{1}).range)),'yticklabel',[pop{d}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{d}.(fields{1}).counts,1)
        text(20,k,num2str(pop{d}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{d}.(fields{1}).range) 0],'w:')
    %       axis('square')
    %             xlabel('time from touch onset (ms)')
    %             title(['theta'])
    
end
