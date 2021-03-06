%% Figure 5 Tuning of neuronal responses to touch order, pretouch velocity
% and maximum curvature. It requires that variables created with the 
% S1datasetWrapper.m and intermediateAnalyasisStructures.mat are in the 
% main workspace.
% Tested on Matlab 2015b
% Samuel Andrew Hires 2016-07-14


%% Figure 5A
% This cell will plot all L4 touch tuning in the style of Figure 5A
% The example figure in paper uses rec = 23.  Other cells can be plotted by
% changing rec to equal other cell numbers.  Code for phase modulated touch
% and printing is present but commented out.  Figure 5B requires that this
% the intermediate variables calculated in this 5A cell are present for rec
% 1:31

 %printdir = 'Z:\Users\Jon\Projects\Characterization\L5b\Figures\LocationTuning\'
for i=1:length(U)
resL45(i).idxT=50;
U{i}.onsetLatency=6;
end
L=U;

if strcmp(U{1}.meta.layer,'L5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I23:J50');
    disp(U{1}.meta.layer);
    layer = 'L5b';
elseif strcmp(U{1}.meta.layer,'L3')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B25:C44');
    disp(U{1}.meta.layer);
    layer = 'L3';
elseif strcmp(U{1}.meta.layer,'L4')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','O23:P28');
    disp(U{1}.meta.layer);
    layer = 'L4';
elseif strcmp(U{1}.meta.layer,'L3Out')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B77:C82');
    disp(U{1}.meta.layer);
    layer = 'L3Out';
else strcmp(U{1}.meta.layer,'L5bOut')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I75:J80');
    disp(U{1}.meta.layer);
    layer = 'L5bOut';
end
txt=txt(~isnan(ndata));
ndata=ndata(~isnan(ndata));

 
for p=1:length(ndata)
    cellcode = txt{p};
    rec = p;
    
    
    % Adaptation across touches
    % Adaptation by ITI
    % Phase at touch
    % Velocity at touch
    % Peak force at touch
    % Direction of touch
        
    trange = [1:resL45(rec).idxT];

    v(rec).spikes = squeeze(L{rec}.R_ntk);
    v(rec).allTouchIdx = find(nansum([L{rec}.S_ctk(9,:,:);L{rec}.S_ctk(12,:,:)]));
        v(rec).allTouchIdx = v(rec).allTouchIdx(v(rec).allTouchIdx<(numel(v(rec).spikes)-151)); %cuts off touches within last 50ms of last trial 
    v(rec).firstTouchIdx = find(L{rec}.S_ctk(9,:,:)==1);
    v(rec).lateTouchIdx = find(L{rec}.S_ctk(12,:,:)==1);
        v(rec).lateTouchIdx = v(rec).lateTouchIdx(v(rec).lateTouchIdx<(numel(v(rec).spikes)-151));%cuts off touches within last 50ms of last trial 
    v(rec).allTouchITI = diff([0; v(rec).allTouchIdx]);
    
    theta = squeeze(L{rec}.S_ctk(1,:,:));
    vel = [zeros(1,L{rec}.k); squeeze(diff( L{rec}.S_ctk(1,:,:)))];
    amp = squeeze(L{rec}.S_ctk(3,:,:));
    setp = squeeze(L{rec}.S_ctk(4,:,:));
    phase = squeeze(L{rec}.S_ctk(5,:,:));
    dk = squeeze(L{rec}.S_ctk(6,:,:));
    M0Adj = squeeze(L{rec}.S_ctk(7,:,:));
    Fax = squeeze(L{rec}.S_ctk(8,:,:));
    tspikesIdx = repmat(v(rec).allTouchIdx,1,length(trange))+repmat(trange+L{rec}.onsetLatency-1,length(v(rec).allTouchIdx),1);
    
    if isnan(tspikesIdx)
        v(rec).sc = 0;
        v(rec).lsc = 0;
    else
        v(rec).sc = sum(v(rec).spikes(tspikesIdx),2);  
        lspikesIdx = repmat(v(rec).lateTouchIdx,1,length(trange))+repmat(trange+L{rec}.onsetLatency-1,length(v(rec).lateTouchIdx),1);
        v(rec).lsc = sum(v(rec).spikes(lspikesIdx),2);
    end
    
    
    % Touch Theta
    v(rec).allTouchTheta = theta(v(rec).allTouchIdx);     
    v(rec).firstTouchTheta = theta(v(rec).firstTouchIdx);     
    [v(rec).theta.sorted v(rec).theta.sortedBy v(rec).theta.binBounds] = binslin(v(rec).allTouchTheta, v(rec).sc, 'equalN',10); 
    
    v(rec).theta.lh = [];
    v(rec).theta.lci = [];
    
    for num = 1:length(v(rec).theta.sorted)
        [lh, lci] = poissfit(v(rec).theta.sorted{num});
        v(rec).theta.lh(num) = lh;
        v(rec).theta.lci(num,:) = lci;
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
    
    

    %% VELOCITY
    
    figure(52);clf;set(gcf,'paperposition',[0 0 4 4]);subplot(5,4,[5 7])
    plot([51 51],[1 size(vel_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(vel_raster'),151), ceil(find(vel_raster')/151) ,'k.','markersize',4)
    ylabel('Pre Velocity')
    axis([1 151 1 size(vel_raster,1)]);
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[],'ytick',[])
    ax1 = gca;
    ax1;hold on
    ax1_pos = get(ax1,'Position');
    axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'yaxislocation','right','Color','none');hold on
    plot(mat2gray(sort_vel(:,2)),1:length(sort_vel),'r')
    axis tight
    norm_vel_range = mat2gray([sort_vel(1,2) 0 sort_vel(end,2)]);
    plot([norm_vel_range(2)]*[1;1],[1 size(sort_vel,1)],'--','color',[1 .5 .5])
    set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_vel(1,2) 0 sort_vel(end,2)]),'xticklabel',round([sort_vel(1,2) 0 (sort_vel(end,2))]))
    %xlabel('Pretouch velocity','color','r')
    
    subplot(5,4,8)
    plot(v(rec).velocity.lh, 1:10,'k')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).velocity.lci,1:10, 'color',[.5 .5 .5])
    
    set(gca,'color','none')
    %    print(gcf,'-depsc2',[printdir 'Raster_Velocity_Cell_' num2str(rec)])
    
    %% TRIAL RASTER
        subplot(5,4,[1 4]);hold on
    title('At touch')
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end
        
    end
    set(gca,'ylim',[0 U{rec}.k]+1)
    %% ADAPTATION
    [sort_adapt,idx_adapt] = sortrows ([v(rec).allTouchIdx, v(rec).touchNumber'],2);
    adapt_raster = zeros(length(idx_adapt),151);
    adapt_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    
    figure(52);subplot(5,4,[9 11])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    ylabel('Touch #')
    plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[],'ytick',[])
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
    %xlabel('Touch order in trial','color','r')
  
    subplot(5,4,12)
    plot(v(rec).adaptation.lh, 1:size(v(rec).adaptation.lh,2),'k-')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[5 10.5],'ytick',[])
    plot(v(rec).adaptation.lci,1:size(v(rec).adaptation.lh,2), 'color',[.5 .5 .5])
    axes(ax1)
    set(gca,'color','none')
    %  print(gcf,'-depsc2',[printdir 'Raster_Adaptation_Cell_' num2str(rec)])
    
%% THETA
    [sort_theta,idx_theta] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchTheta],2);
    [sort_theta_1st,idx_theta_1st] = sortrows ([v(rec).firstTouchIdx, v(rec).firstTouchTheta],2);

    theta_raster = zeros(length(idx_theta),151);
    theta_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_theta),1,151)+repmat([-50:100],length(idx_theta),1));
    
    theta_raster_1st = v(rec).spikes(repmat(v(rec).firstTouchIdx(idx_theta_1st),1,151)+repmat([-50:100],length(idx_theta_1st),1));
    
    
    figure(52);subplot(5,4,[13 15])
    plot([51 51],[1 size(theta_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(theta_raster'),151), ceil(find(theta_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(theta_raster,1)]);
    colormap([1 1 1;0 0 0])
    ylabel('Theta')
    box off
    set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[],'ytick',[])
    ax1 = gca;
    ax1;hold on
    ax1_pos = get(ax1,'Position');
    axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'yaxislocation','right','Color','none');hold on
    plot(mat2gray(sort_theta(:,2)),1:length(sort_theta),'r')
    axis tight
    norm_adapt_range = mat2gray([sort_theta(1,2) sort_theta(end,2)]);
    plot([norm_adapt_range(2)]*[1;1],[1 size(sort_theta,1)],'--','color',[1 .5 .5])
    set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_theta(1,2) sort_theta(end,2)]),'xticklabel',round([sort_theta(1,2) (sort_theta(end,2))]))
  
    figure(52);subplot(5,4,16)
    plot(v(rec).theta.lh, 1:size(v(rec).theta.lh,2),'k-')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).theta.lci,1:size(v(rec).theta.lh,2), 'color',[.5 .5 .5])
    axes(ax1)
    set(gca,'color','none')

%% KAPPA
    
    k_range = [1:resL45(rec).idxT]; % Pre touch indicies
    k_Idx = repmat(v(rec).allTouchIdx,1,length(k_range))+repmat(k_range,length(v(rec).allTouchIdx),1);
    allTouchK = dk(k_Idx)-repmat(dk(v(rec).allTouchIdx+1),1,length(k_range));
    [~, kmax_idx] = max(abs(allTouchK'));
    v(rec).allTouchK = [];
    for i = 1:length(kmax_idx)
        v(rec).allTouchK(i) = allTouchK(i,kmax_idx(i));
    end
    v(rec).allTouchK = v(rec).allTouchK';
    [sort_k,idx_k] = sortrows([v(rec).allTouchIdx, v(rec).allTouchK],2);
    k_raster = zeros(length(idx_k),151);
    k_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_k),1,151)+repmat([-50:100],length(idx_k),1));
    
    figure(52);subplot(5,4,[1 3]) 
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    ylabel('Kappa')
    plot(mod(find(k_raster'),151), ceil(find(k_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 101],'xtick',[],'xticklabel',[],'ytick',[])
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
    
    [v(rec).k.sorted v(rec).k.sortedBy v(rec).k.binBounds]=binslin(v(rec).allTouchK, v(rec).sc, 'equalN',10);
    v(rec).k.means = cellfun(@nanmean, v(rec).k.sorted);
    v(rec).k.lh = [];
    v(rec).k.lci = [];
    for num = 1:length(v(rec).k.sorted)
        [lh, lci] = poissfit(v(rec).k.sorted{num});
        v(rec).k.lh(num) = lh;
        v(rec).k.lci(num,:) = lci;
    end
    
    subplot(5,4,[4])
    plot(v(rec).k.lh, 1:10,'k')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).k.lci,1:10, 'color',[.5 .5 .5])
    set(gca,'color','none') 
        %print(gcf,'-depsc2',[printdir 'Raster_Combo_Cell_' num2str(rec)])    
%     subplot(5,4,[1 3]) 
%      plot([51 51],[1 size(theta_raster_1st,1)],'--','color',[.5 .5 .5])
%     hold on
%     plot(mod(find(theta_raster_1st'),151), ceil(find(theta_raster_1st')/151) ,'k.','markersize',4)
%     axis([1 151 1 size(theta_raster_1st,1)]);
%     colormap([1 1 1;0 0 0])
%     box off
%     set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
%     ax1 = gca;
%     ax1;hold on
%     ax1_pos = get(ax1,'Position');
%     axes('Position',ax1_pos,...
%         'XAxisLocation','top',...
%         'yaxislocation','right','Color','none');hold on
%     plot(mat2gray(sort_theta_1st(:,2)),1:length(sort_theta_1st),'r')
%     axis tight
%     norm_adapt_range = mat2gray([sort_theta_1st(1,2) sort_theta_1st(end,2)]);
%     plot([norm_adapt_range(2)]*[1;1],[1 size(sort_theta_1st,1)],'--','color',[1 .5 .5])
%     set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_theta_1st(1,2) sort_theta_1st(end,2)]),'xticklabel',round([sort_theta_1st(1,2) (sort_theta_1st(end,2))]))

    %print(gcf,'-depsc2',[printdir 'Raster_Adaptation_Cell_' num2str(rec)])
      % print(gcf,'-depsc2',[printdir 'Location_Tuning_JC_L3Cell_' num2str(rec)])
 %%
%  %%  Conditioned
%     conIdx = find(v(rec).allTouchVel >0);
%     v(rec).conTouchIdx = v(rec).allTouchIdx(conIdx);
%     v(rec).conTouchTheta = v(rec).allTouchTheta(conIdx);
%      
%     [sort_theta_con,idx_theta_con] = sortrows ([v(rec).conTouchIdx, v(rec).conTouchTheta],2);
%  
%  
%         theta_raster_con = v(rec).spikes(repmat(v(rec).conTouchIdx(idx_theta_con),1,151)+repmat([-50:100],length(idx_theta_con),1));
% 
%     figure(56);clf;
%     plot([51 51],[1 size(theta_raster_con,1)],'--','color',[.5 .5 .5])
%     hold on
%     plot(mod(find(theta_raster_con'),151), ceil(find(theta_raster_con')/151) ,'k.','markersize',4)
%     axis([1 151 1 size(theta_raster_con,1)]);
%     colormap([1 1 1;0 0 0])
%     box off
%     set(gca,'xlim',[26 101],'xtick',[26 51 101],'xticklabel',[-25 0 50],'ytick',[])
%     ax1 = gca;
%     ax1;hold on
%     ax1_pos = get(ax1,'Position');
%     axes('Position',ax1_pos,...
%         'XAxisLocation','top',...
%         'yaxislocation','right','Color','none');hold on
%     plot(mat2gray(sort_theta_con(:,2)),1:length(sort_theta_con),'r')
%     axis tight
%     norm_adapt_range = mat2gray([sort_theta_con(1,2) sort_theta_con(end,2)]);
%     plot([norm_adapt_range(2)]*[1;1],[1 size(sort_theta_con,1)],'--','color',[1 .5 .5])
%     set(gca,'ytick',[],'xcolor','r','xtick',mat2gray([sort_theta_con(1,2) sort_theta_con(end,2)]),'xticklabel',round([sort_theta_con(1,2) (sort_theta_con(end,2))]))
%     xlabel('Theta at touch onset','color','r')

 
 %% PHASE  Optional section
        v(rec).allTouchPhase = phase(v(rec).allTouchIdx)
        [sort_phase,idx_phase] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchPhase],2);
        phase_raster = zeros(length(idx_phase),151);
        phase_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_phase),1,151)+repmat([-50:100],length(idx_phase),1));
    
        figure(52);subplot(5,4,[17 19])
        plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
        hold on
        plot(mod(find(phase_raster'),151), ceil(find(phase_raster')/151) ,'k.','markersize',4)
        axis([1 151 1 size(adapt_raster,1)]);
        xlabel('Time from touch (ms)')
        ylabel('Phase')
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
        %xlabel('Phase at touch','color','r')
        figure(52);subplot(5,4,20)
        plot(v(rec).phase.lh, 1:12,'k')
        hold on
        box off
        set(gca,'ylim',[0.5 12.5],'ytick',[])
        xlabel('spks/touch')
        plot(v(rec).phase.lci,1:12, 'color',[.5 .5 .5])
        axes(ax1)
        set(gca,'color','none')
     %print(gcf,'-depsc2',[printdir 'Raster_Phase_Cell_' num2str(rec)])
    print(gcf,'-dpng',['Z:\Users\Jon\Projects\Characterization\' layer '\Figures\' cellcode '_' 'Variables_at_touch_Sorted_Raster'])
 
    %%  Modulation depth
    
    %%w(rec).vel    = (max(v(rec).velocity.lh)-min(v(rec).velocity.lh))/(max(v(rec).velocity.lh)+min(v(rec).velocity.lh));
    %w(rec).adapt  = (max(v(rec).adaptation.lh(1:10))-min(v(rec).adaptation.lh(1:10)))/(max(v(rec).adaptation.lh(1:10))+min(v(rec).adaptation.lh(1:10)));
    %w(rec).M0     = (max(v(rec).k.lh)-min(v(rec).k.lh))/(max(v(rec).k.lh)+min(v(rec).k.lh));
     %w(rec).phase  = (max(v(rec).phase.lh)-min(v(rec).phase.lh))/(max(v(rec).phase.lh)+min(v(rec).phase.lh));
    
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
