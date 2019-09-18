%% Figure 5 Tuning of neuronal responses to touch order, pretouch velocity
% and maximum curvature. It requires that variables created with the 
% S1datasetWrapper.m and intermediateAnalyasisStructures.mat are in the 
% main workspace.
% Tested on Matlab 2015b
% Samuel Andrew Hires 2016-07-14


%% Figure 5A

selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));

for rec = 1:length(selectedCells)
    array=U{selectedCells(rec)}
    array.onsetLatency=array.meta.touchProperties.responseWindow(1);
    trange=1:diff(array.meta.touchProperties.responseWindow);
    % Adaptation across touches
    % Adaptation by ITI
    % Phase at touch
    % Velocity at touch
    % Peak force at touch
    % Direction of touch
    
    v(rec).spikes = squeeze(array.R_ntk);
    v(rec).allTouchIdx = find(nansum([array.S_ctk(9,:,:);array.S_ctk(12,:,:)]));
    v(rec).firstTouchIdx = find(array.S_ctk(9,:,:)==1);
    v(rec).lateTouchIdx = find(array.S_ctk(12,:,:)==1);
    v(rec).allTouchITI = diff([0; v(rec).allTouchIdx]);
    
    vel = [zeros(1,array.k); squeeze(diff( array.S_ctk(1,:,:)))];
    theta=squeeze(array.S_ctk(1,:,:));
    amp = squeeze(array.S_ctk(3,:,:));
    setp = squeeze(array.S_ctk(4,:,:));
    phase = squeeze(array.S_ctk(5,:,:));
    dk = squeeze(array.S_ctk(6,:,:));
    M0Adj = squeeze(array.S_ctk(7,:,:));
    Fax = squeeze(array.S_ctk(8,:,:));
    
    %trange = [1:resL45(rec).idxT];
    tspikesIdx = repmat(v(rec).allTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(v(rec).allTouchIdx),1); %index for touch + trange
    
    if isnan(tspikesIdx)
        v(rec).sc = 0;
        v(rec).lsc = 0;
    else
        v(rec).sc = sum(v(rec).spikes(tspikesIdx),2);  
        lspikesIdx = repmat(v(rec).lateTouchIdx,1,length(trange))+repmat(trange+array.onsetLatency-1,length(v(rec).lateTouchIdx),1);
        v(rec).lsc = sum(v(rec).spikes(lspikesIdx),2);
    end
    
    % Touch Adaptation
    v(rec).touchNumber = [];
    for i = 1:length(v(rec).firstTouchIdx)-1;
        v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(i) & v(rec).allTouchIdx < v(rec).firstTouchIdx(i+1)));
    end
    v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(end)));
    v(rec).lateTouchNumber = v(rec).touchNumber(v(rec).touchNumber>1);
    v(rec).adaptation.lh = [];
    v(rec).adaptation.lci = [];
 
    [lh, lci] = poissfit(v(rec).sc(v(rec).touchNumber == 1));
    v(rec).adaptation.lh(1) = lh;
    v(rec).adaptation.lci(1,:) = lci;
    [v(rec).adaptation.sorted v(rec).adaptation.sortedBy v(rec).adaptation.binBounds]=binslin(v(rec).lateTouchNumber, v(rec).lsc, 'equalN',9);
    
    for num = 1:length(v(rec).adaptation.sorted)
        [lh, lci] = poissfit(v(rec).adaptation.sorted{num});
        v(rec).adaptation.lh(num+1) = lh;
        v(rec).adaptation.lci(num+1,:) = lci;
    end
    
    %Touch Velocity
    vel_range = [-4:-1]; % Pre touch indicies
    vel_Idx = repmat(v(rec).allTouchIdx,1,length(vel_range))+repmat(vel_range,length(v(rec).allTouchIdx),1);
    v(rec).allTouchVel = nanmean(vel(vel_Idx),2);
    [v(rec).velocity.sorted v(rec).velocity.sortedBy v(rec).velocity.binBounds]=binslin(v(rec).allTouchVel, v(rec).sc, 'equalN',10);
    v(rec).velocity.means = cellfun(@nanmean, v(rec).velocity.sorted);
    v(rec).velocity.lh = [];
    v(rec).velocity.lci = [];
 
    for num = 1:length(v(rec).velocity.sorted)
        [lh, lci] = poissfit(v(rec).velocity.sorted{num});
        v(rec).velocity.lh(num) = lh;
        v(rec).velocity.lci(num,:) = lci;
    end
    
    [sort_vel,idx_vel] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchVel],2);
    vel_raster = zeros(length(idx_vel),151);
    vel_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_vel),1,151)+repmat([-50:100],length(idx_vel),1));
    
    % Touch Theta
    v(rec).allTouchTheta = theta(v(rec).allTouchIdx); %find all theta at touches
   [v(rec).theta.sorted, v(rec).theta.sortedBy ,v(rec).theta.binBounds]=binslin(v(rec).allTouchTheta, v(rec).sc, 'equalE',41,-50,50);
   v(rec).theta.means = cellfun(@nanmean, v(rec).theta.sorted); %mean spike count in each bin 
     v(rec).theta.lh = [];
    v(rec).theta.lci = [];
    
    for num = 1:length(v(rec).theta.sorted) %poissfit to identify likelihood of event occuring
        [lh, lci] = poissfit(v(rec).theta.sorted{num});
        v(rec).theta.lh(num) = lh; %lambdahat = MLE mean for poissfit
        v(rec).theta.lci(num,:) = lci; %lambda confiddence interval for lambda mean
    end
    
   [sort_theta,idx_theta] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchTheta],2);
    
   theta_raster = zeros(length(idx_theta),151); %this is for building raster of theta and the spikes -50ms:100ms post touch
   theta_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_theta),1,151)+repmat([-50:100],length(idx_theta),1));
    
    %% ADAPTATION
    [sort_adapt,idx_adapt] = sortrows ([v(rec).allTouchIdx, v(rec).touchNumber'],2);
    adapt_raster = zeros(length(idx_adapt),151);
    adapt_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    
    figure(52);subplot(3,4,[1 3])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
    ax1 = gca;
    ax1;hold on
    ax1_pos = get(ax1,'Position');
    axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'yaxislocation','right','Color','none');hold on
    plot(mat2gray(sort_adapt(:,2)),1:length(sort_adapt),'r')
    axis tight
    norm_adapt_range = mat2gray([sort_adapt(1,2) sort_adapt(end,2)]);
    plot([norm_adapt_range(2)]*[1;1],[1 size(sort_adapt,1)],'--','color',[1 .5 .5])
    set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_adapt(1,2) sort_adapt(end,2)]),'xticklabel',round([sort_adapt(1,2) (sort_adapt(end,2))]))
    xlabel('Touch order in trial','color','r')
  
    subplot(3,4,4)
    plot(v(rec).adaptation.lh, 1:size(v(rec).adaptation.lh,2),'k-')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).adaptation.lci,1:size(v(rec).adaptation.lh,2), 'color',[.5 .5 .5])
    axes(ax1)
    set(gca,'color','none')
    %  print(gcf,'-depsc2',[printdir 'Raster_Adaptation_Cell_' num2str(rec)])
    
    %% KAPPA
    
    %k_range = [1:resL45(rec).idxT]; % Pre touch indicies
    k_Idx = repmat(v(rec).allTouchIdx,1,length(k_range))+repmat(k_range,length(v(rec).allTouchIdx),1);
    allTouchK = dk(k_Idx)-repmat(dk(v(rec).allTouchIdx+1),1,length(k_range));
    [~, kmax_idx] = max(abs(allTouchK'));
    v(rec).allTouchK = [];
    for i = 1:length(kmax_idx)
        v(rec).allTouchK(i) = allTouchK(i,kmax_idx(i));
    end
    v(rec).allTouchK = v(rec).allTouchK'
    [sort_k,idx_k] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchK],2);
    k_raster = zeros(length(idx_k),151);
    k_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_k),1,151)+repmat([-50:100],length(idx_k),1));
    
    figure(52);subplot(3,4,[9 11])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(k_raster'),151), ceil(find(k_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    xlabel('Time from touch (ms)')
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
    ax1 = gca;
    ax1;hold on
    ax1_pos = get(ax1,'Position');
    axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'yaxislocation','right','Color','none');hold on
    plot(sort_k(:,2),1:length(sort_k),'r')
    axis tight
    norm_k_range = mat2gray([sort_k(1,2) 0 sort_k(end,2)]);
    plot([0 0],[1 size(sort_k,1)],'--','color',[1 .5 .5])
    set(gca,'ytick',[],'xcolor','r','xtick',[min(sort_k(:,2)) 0 max(sort_k(:,2))],'xticklabel',[sort_k(1,2) 0 (sort_k(end,2))])
    xlabel('Max Touch Kappa')
    
    [v(rec).k.sorted v(rec).k.sortedBy v(rec).k.binBounds]=binslin(v(rec).allTouchK, v(rec).sc, 'equalN',10);
    v(rec).k.means = cellfun(@nanmean, v(rec).k.sorted);
    v(rec).k.lh = [];
    v(rec).k.lci = [];
    for num = 1:length(v(rec).k.sorted)
        [lh, lci] = poissfit(v(rec).k.sorted{num});
        v(rec).k.lh(num) = lh;
        v(rec).k.lci(num,:) = lci;
    end
    
    subplot(3,4,12)
    plot(v(rec).k.lh, 1:10,'k')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).k.lci,1:10, 'color',[.5 .5 .5])
    xlabel('spks/touch')
    set(gca,'color','none')
    %     print(gcf,'-depsc2',[printdir 'Raster_Combo_Cell_' num2str(rec)])
    
    %% PHASE  Optional section
        v(rec).allTouchPhase = phase(v(rec).allTouchIdx)
        [sort_phase,idx_phase] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchPhase],2);
        phase_raster = zeros(length(idx_phase),151);
        phase_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_phase),1,151)+repmat([-50:100],length(idx_phase),1));
    
        figure(52);subplot(4,4,[9 11])
        plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
        hold on
        plot(mod(find(phase_raster'),151), ceil(find(phase_raster')/151) ,'k.','markersize',4)
        axis([1 151 1 size(adapt_raster,1)]);
        xlabel('Time from touch (ms)')
        ylabel('Touches sorted by phase')
        colormap([1 1 1;0 0 0])
        box off
        set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
        ax1 = gca;
        ax1;hold on
        ax1_pos = get(ax1,'Position');
        axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'yaxislocation','right','Color','none');hold on
        plot(sort_phase(:,2),1:length(sort_phase),'r')
        axis tight
        set(gca,'xlim',[-pi pi],'xcolor','r','xtick',[-pi 0 pi],'xticklabel',{'-pi', 0, 'pi'},'ytick',[])
        norm_phase_range = mat2gray([sort_phase(1,2) 0 sort_phase(end,2)]);
        plot([0 0],[1 size(sort_phase,1)],'--','color',[1 .5 .5])
        [v(rec).phase.sorted v(rec).phase.sortedBy v(rec).phase.binBounds]=binslin(v(rec).allTouchPhase, v(rec).sc, 'equalN',12)
        v(rec).phase.means = cellfun(@nanmean, v(rec).phase.sorted);
    
        v(rec).phase.lh = [];
        v(rec).phase.lci = [];
        for num = 1:length(v(rec).phase.sorted)
            [lh, lci] = poissfit(v(rec).phase.sorted{num})
            v(rec).phase.lh(num) = lh;
            v(rec).phase.lci(num,:) = lci;
        end
        xlabel('Phase at touch','color','r')
        subplot(4,4,12)
        plot(v(rec).phase.lh, 1:12,'k')
        hold on
        box off
        set(gca,'ylim',[0.5 12.5],'ytick',[])
        plot(v(rec).phase.lci,1:12, 'color',[.5 .5 .5])
        xlabel('spks/touch')
        axes(ax1)
        set(gca,'color','none')
    %  print(gcf,'-depsc2',[printdir 'Raster_Phase_Cell_' num2str(rec)])
    
    %%  Modulation depth
    
    w(rec).vel    = (max(v(rec).velocity.lh)-min(v(rec).velocity.lh))/(max(v(rec).velocity.lh)+min(v(rec).velocity.lh));
    w(rec).adapt  = (max(v(rec).adaptation.lh(1:10))-min(v(rec).adaptation.lh(1:10)))/(max(v(rec).adaptation.lh(1:10))+min(v(rec).adaptation.lh(1:10)));
    w(rec).M0     = (max(v(rec).k.lh)-min(v(rec).k.lh))/(max(v(rec).k.lh)+min(v(rec).k.lh));
%     w(rec).phase  = (max(v(rec).phase.lh)-min(v(rec).phase.lh))/(max(v(rec).phase.lh)+min(v(rec).phase.lh));
    pause
end


%% Figure 5B % Population data
% The sorted order of cells is slightly different from the original publication 
% due to a refinement of the set of performing trials in some sessions

vel_map = [];
adp_map = [];
kap_map = [];

for rec = 1:31
    vel_map(rec,:) = v(rec).velocity.lh/max(v(rec).velocity.lh);
    adp_map(rec,:) = v(rec).adaptation.lh/max(v(rec).adaptation.lh);
    kap_map(rec,:) = v(rec).k.lh/max(v(rec).k.lh);
end

figure(10);
clf;
set(gcf,'paperposition',[0 0 1.5 4])
[sortMap sortMapIdx] = sort(kap_map*linspace(-1,1,10)')
subplot(3,1,1)
image(adp_map(sortMapIdx,:)'*256)
colormap(gray(256))
axis xy
set(gca,'xtick',[],'ytick',[])
ylabel('Touch number')

subplot(3,1,2)
image(vel_map(sortMapIdx,:)'*256)
colormap(gray(256))
axis xy
set(gca,'xtick',[],'ytick',[])
ylabel('Velocity')

subplot(3,1,3)
image(kap_map(sortMapIdx,:)'*256)
colormap(gray(256))
xlabel('Cells')
axis xy
set(gca,'xtick',[],'ytick',[])
ylabel('Kappa')

%print(gcf,'-depsc2',[printdir 'PopulationVarianceByPredictor'])

%% Figure 5C
figure(11);clf
set(gcf,'paperposition',[0 0 2.5 1.5])

image([w(sortMapIdx).adapt;w(sortMapIdx).vel;w(sortMapIdx).M0]*256)
colormap(gray(256));
title('Modulation depth')
set(gca,'xtick',[],'ytick',[])
h_cb = colorbar;
set(h_cb,'ytick',[1 128 256],'yticklabel',[0 .5 1])

%print(gcf,'-depsc2',[printdir 'PopulationModulationDepth'])

