printdir = 'Z:\Users\Jon\DATA\Figures\work\'
for rec = 1
    
    % Adaptation across touches
    % Adaptation by ITI
    % Phase at touch
    % Velocity at touch
    % Peak force at touch
    % Direction of touch
    
    v(rec).spikes = squeeze(U{rec}.R_ntk);
    
    v(rec).allTouchIdx = find(nansum([U{rec}.S_ctk(9,:,:);U{rec}.S_ctk(12,:,:)]));
    v(rec).firstTouchIdx = find(U{rec}.S_ctk(9,:,:)==1);
    v(rec).lateTouchIdx = find(U{rec}.S_ctk(12,:,:)==1);
    v(rec).allTouchITI = diff([0; v(rec).allTouchIdx]);
    
    vel = squeeze(U{rec}.S_ctk(2,:,:));
    amp = squeeze(U{rec}.S_ctk(3,:,:));
    setp = squeeze(U{rec}.S_ctk(4,:,:));
    phase = squeeze(U{rec}.S_ctk(5,:,:));
    dk = squeeze(U{rec}.S_ctk(6,:,:));
    M0Adj = squeeze(U{rec}.S_ctk(7,:,:));
    Fax = squeeze(U{rec}.S_ctk(8,:,:));
    
    %trange = [1:resL45(rec).idxT];
    trange = [34:60];
    tspikesIdx = repmat(v(rec).allTouchIdx,1,length(trange))+repmat(trange+U{rec}.meta.onsetLatency-1,length(v(rec).allTouchIdx),1);
    if isnan(tspikesIdx)
        v(rec).sc = 0;
        v(rec).lsc = 0;
    else
        
    v(rec).sc = sum(v(rec).spikes(tspikesIdx),2)
    
    lspikesIdx = repmat(v(rec).lateTouchIdx,1,length(trange))+repmat(trange+U{rec}.meta.onsetLatency-1,length(v(rec).lateTouchIdx),1);
    v(rec).lsc = sum(v(rec).spikes(lspikesIdx),2)
    end
    % Touch Adaptation
    v(rec).touchNumber = []
    for i = 1:length(v(rec).firstTouchIdx)-1;
        v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(i) & v(rec).allTouchIdx < v(rec).firstTouchIdx(i+1)))
    end
    v(rec).touchNumber = cat(2,v(rec).touchNumber,1:sum(v(rec).allTouchIdx >= v(rec).firstTouchIdx(end)));
    v(rec).lateTouchNumber = v(rec).touchNumber(v(rec).touchNumber>1);

    v(rec).adaptation.lh = [];
    v(rec).adaptation.lci = [];
    
    [lh, lci] = poissfit(v(rec).sc(v(rec).touchNumber == 1))
    v(rec).adaptation.lh(1) = lh;
    v(rec).adaptation.lci(1,:) = lci;
    
    [v(rec).adaptation.sorted v(rec).adaptation.sortedBy v(rec).adaptation.binBounds]=binslin(v(rec).lateTouchNumber, v(rec).lsc, 'equalN',9)

    for num = 1:length(v(rec).adaptation.sorted)
        [lh, lci] = poissfit(v(rec).adaptation.sorted{num})
        v(rec).adaptation.lh(num+1) = lh;
        v(rec).adaptation.lci(num+1,:) = lci;
    end
    
    %Touch Velocity
    vel_range = [-4:-1]; % Pre touch indicies
    vel_Idx = repmat(v(rec).allTouchIdx,1,length(vel_range))+repmat(vel_range,length(v(rec).allTouchIdx),1);
    v(rec).allTouchVel = nanmean(vel(vel_Idx),2)
    
    
    [v(rec).velocity.sorted v(rec).velocity.sortedBy v(rec).velocity.binBounds]=binslin(v(rec).allTouchVel, v(rec).sc, 'equalN',10)
    v(rec).velocity.means = cellfun(@nanmean, v(rec).velocity.sorted);
    
    v(rec).velocity.lh = [];
    v(rec).velocity.lci = [];
    for num = 1:length(v(rec).velocity.sorted)
        [lh, lci] = poissfit(v(rec).velocity.sorted{num})
        v(rec).velocity.lh(num) = lh;
        v(rec).velocity.lci(num,:) = lci;
    end
    [sort_vel,idx_vel] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchVel],2);
    vel_raster = zeros(length(idx_vel),151);
    vel_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_vel),1,151)+repmat([-50:100],length(idx_vel),1));
    
    
    h_1 = figure(1);clf;subplot(3,2,1)
    plot(v(rec).adaptation.lci,'color',[.5 .5 1])
    hold on
    plot(v(rec).adaptation.lh)
    title('Adaptation')
    
    subplot(3,2,2)
    plot((v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2,v(rec).velocity.lh, 'color','b')
    hold on
    plot(repmat((v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2,2,1)',v(rec).velocity.lci, 'color',[.5 .5 1])
    title('Velocity')
    
    %set(ax2,'xtick',mat2gray([sort_vel(1,2) 0 sort_vel(end,2)]),'xticklabel',[sort_vel(1,2) 0 sort_vel(end,2)])
    
    %  Raster plots
    %% VELOCITY
    figure(2);clf;set(gcf,'paperposition',[0 0 4 4]);subplot(3,4,[5 7])
    %imagesc(vel_raster)
    plot([51 51],[1 size(vel_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(vel_raster'),151), ceil(find(vel_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(vel_raster,1)]);
    %xlabel('Time from touch (ms)')
    %ylabel('Sorted touches')
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 126],'xtick',[26 51 76 100 126],'xticklabel',[-25 0 25 50 75 100],'ytick',[])
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
    xlabel('Pretouch velocity','color','r')
    subplot(3,4,8)
    plot(v(rec).velocity.lh, 1:10,'k')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).velocity.lci,1:10, 'color',[.5 .5 .5])
    %xlabel('spks/touch')
    %    axes(ax1)
    set(gca,'color','none')
     print(gcf,'-depsc2',[printdir 'Raster_Velocity_Cell_' num2str(rec)])
    
    %% ADAPTATION
    [sort_adapt,idx_adapt] = sortrows ([v(rec).allTouchIdx, v(rec).touchNumber'],2);
    adapt_raster = zeros(length(idx_adapt),151);
    adapt_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_adapt),1,151)+repmat([-50:100],length(idx_adapt),1));
    
    figure(2);subplot(3,4,[1 3])
    %imagesc(vel_raster)
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(adapt_raster'),151), ceil(find(adapt_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    %xlabel('Time from touch (ms)')
    %ylabel('Touches sorted by order')
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 126],'xtick',[26 51 76 100 126],'xticklabel',[-25 0 25 50 75 100],'ytick',[])
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
    %xlabel('spks/touch')
    axes(ax1)
    set(gca,'color','none')
      print(gcf,'-depsc2',[printdir 'Raster_Adaptation_Cell_' num2str(rec)])
    
    
    
    %% PHASE
    v(rec).allTouchPhase = phase(v(rec).allTouchIdx)
    [sort_phase,idx_phase] = sortrows ([v(rec).allTouchIdx, v(rec).allTouchPhase],2);
    phase_raster = zeros(length(idx_phase),151);
    phase_raster = v(rec).spikes(repmat(v(rec).allTouchIdx(idx_phase),1,151)+repmat([-50:100],length(idx_phase),1));
    
    figure(1);subplot(4,4,[9 11])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(phase_raster'),151), ceil(find(phase_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    xlabel('Time from touch (ms)')
    ylabel('Touches sorted by phase')
    colormap([1 1 1;0 0 0])
    box off
    set(gca,'xlim',[26 126],'xtick',[26 51 76 100 126],'xticklabel',[-25 0 25 50 75 100],'ytick',[])
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
    plot(v(rec).phase.lh, 1:12,'k')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 12.5],'ytick',[])
    plot(v(rec).phase.lci,1:12, 'color',[.5 .5 .5])
    xlabel('spks/touch')
    axes(ax1)
    set(gca,'color','none')
      print(gcf,'-depsc2',[printdir 'Raster_Phase_Cell_' num2str(rec)])
    
    
    %% KAPPA
    
    k_range = 1:max(trange); % Pre touch indicies
    k_Idx = repmat(v(rec).allTouchIdx,1,length(k_range))+repmat(k_range,length(v(rec).allTouchIdx),1);
    %k_Idx = k_Idx(1:end-1,:);
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
    
    figure(2);subplot(3,4,[9 11])
    plot([51 51],[1 size(adapt_raster,1)],'--','color',[.5 .5 .5])
    hold on
    plot(mod(find(k_raster'),151), ceil(find(k_raster')/151) ,'k.','markersize',4)
    axis([1 151 1 size(adapt_raster,1)]);
    xlabel('Time from touch (ms)')
    %ylabel('Touches sorted by M0')
    colormap([1 1 1;0 0 0])
    
    
    box off
    set(gca,'xlim',[26 126],'xtick',[26 51 76 100 126],'xticklabel',[-25 0 25 50 75 100],'ytick',[])
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
    set(gca,'ytick',[],'xcolor','r','xtick',[sort_k(1,2) 0 sort_k(end,2)],'xticklabel',[sort_k(1,2) 0 (sort_k(end,2))])
    xlabel('Max Touch Kappa')
    
    [v(rec).k.sorted v(rec).k.sortedBy v(rec).k.binBounds]=binslin(v(rec).allTouchK, v(rec).sc, 'equalN',10)
    v(rec).k.means = cellfun(@nanmean, v(rec).k.sorted);
    v(rec).k.lh = [];
    v(rec).k.lci = [];
    for num = 1:length(v(rec).k.sorted)
        [lh, lci] = poissfit(v(rec).k.sorted{num})
        v(rec).k.lh(num) = lh;
        v(rec).k.lci(num,:) = lci;
    end
    
    subplot(3,4,12)
    plot(v(rec).k.lh, 1:10,'k')%(v(rec).velocity.binBounds(1:end-1)+v(rec).velocity.binBounds(2:end))/2, 'color','b')
    hold on
    box off
    set(gca,'ylim',[0.5 10.5],'ytick',[])
    plot(v(rec).k.lci,1:10, 'color',[.5 .5 .5])
    xlabel('spks/touch')
    %    axes(ax1)
    set(gca,'color','none')
    set(gcf,'paperposition',[0 0 8 6])
    
    print(gcf,'-depsc2',[printdir 'Raster_M0_Cell_' num2str(rec)])
      print(gcf,'-depsc2',[printdir 'Raster_Combo_Cell_' num2str(rec)])
      print(gcf,'-dtiff',[printdir 'Raster_Combo_Cell_' num2str(rec)])
    
    
    
    %%  Modulation depth
    
    w(rec).vel    = (max(v(rec).velocity.lh)-min(v(rec).velocity.lh))/(max(v(rec).velocity.lh)+min(v(rec).velocity.lh));
    w(rec).adapt  = (max(v(rec).adaptation.lh(1:10))-min(v(rec).adaptation.lh(1:10)))/(max(v(rec).adaptation.lh(1:10))+min(v(rec).adaptation.lh(1:10)));
    w(rec).M0     = (max(v(rec).k.lh)-min(v(rec).k.lh))/(max(v(rec).k.lh)+min(v(rec).k.lh));
    w(rec).phase  = (max(v(rec).phase.lh)-min(v(rec).phase.lh))/(max(v(rec).phase.lh)+min(v(rec).phase.lh));
    
end
%%

ww = w([1]);

figure(3);clf;set(gcf,'paperposition',[0 0 6 3])
subplot(1,3,1)
set(gca,'colororder',jet(31))
plot([[ww.vel];zeros(1,31)],[[ww.M0];zeros(1,31)],'o','markersize',4)
hold on
plot([0 1],[0 1],'k:')
xlabel('Velocity')
ylabel('Kappa')
axis square

subplot(1,3,2)
set(gca,'colororder',jet(31))
plot([[ww.M0];zeros(1,31)],[[ww.adapt];zeros(1,31)],'o','markersize',4)
hold on
plot([0 1],[0 1],'k:')
xlabel('Kappa')
ylabel('Adaptation')
axis square

subplot(1,3,3)
set(gca,'colororder',jet(31))
plot([[ww.adapt];zeros(1,31)],[[ww.vel];zeros(1,31)],'o','markersize',4)
hold on
plot([0 1],[0 1],'k:')
xlabel('Adaptation')
ylabel('Velocity')
axis square
    print(gcf,'-depsc2',[printdir 'PopulationVarianceByPredictor' num2str(rec)])


%%
vel_map = [];
adp_map = [];
kap_map = [];
vv = v([1:17 38:51])
for rec = 1:31
    vel_map(rec,:) = vv(rec).velocity.lh/max(vv(rec).velocity.lh);
    adp_map(rec,:) = vv(rec).adaptation.lh/max(vv(rec).adaptation.lh);
    kap_map(rec,:) = vv(rec).k.lh/max(vv(rec).k.lh)
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

print(gcf,'-depsc2',[printdir 'PopulationVarianceByPredictor'])

%%
figure(11);clf
set(gcf,'paperposition',[0 0 2.5 1.5])

image([ww(sortMapIdx).adapt;ww(sortMapIdx).vel;ww(sortMapIdx).M0]*256)
colormap(gray(256))
title('Modulation depth')
set(gca,'xtick',[],'ytick',[])
h_cb = colorbar
set(h_cb,'ytick',[1 128 256],'yticklabel',[0 .5 1])

print(gcf,'-depsc2',[printdir 'PopulationModulationDepth'])

% vel_range = [-4:-1]; % Pre touch indicies
% Fax_range = [-50:150];
% vel_Idx = repmat(v(rec).allTouchIdx,1,length(vel_range))+repmat(vel_range,length(v(rec).allTouchIdx),1);
% Fax_Idx = repmat(v(rec).allTouchIdx,1,length(Fax_range))+repmat(Fax_range,length(v(rec).allTouchIdx),1);
% Fax_Idx_ft = repmat(v(rec).firstTouchIdx,1,length(Fax_range))+repmat(Fax_range,length(v(rec).firstTouchIdx),1);
% v(rec).pcth = nansum(v(rec).spikes(Fax_Idx_ft))
%
% v(rec).allTouchVel = nanmean(vel(vel_Idx),2)
% %v(rec).alltouchPhase = phase(v(rec).allTouchIdx);
% %v(rec).allpeakM0 = M0Adj(tspikesIdx)\
% v(rec).ftFax = Fax(Fax_Idx_ft);
% v(rec).allFax = Fax(Fax_Idx);
%
% ft_idx = find(v(rec).allTouchITI>1000);
% lt_idx = find(v(rec).allTouchITI<1000);
%
%
%
% figure(1);clf;
% subplot(3,3,1);
% [sorted sortedBy binBounds]=binslin(v(rec).alltouchPhase(lt_idx), s.sc(lt_idx), 'equalN',20)
% plot((binBounds(1:end-1)+binBounds(2:end))/2, cellfun(@mean,sorted))
% [sorted sortedBy binBounds]=binslin(v(rec).alltouchPhase(ft_idx), s.sc(ft_idx), 'equalN', 20)
% plot((binBounds(1:end-1)+binBounds(2:end))/2, cellfun(@mean,sorted),'r')
% plot(linspace(-pi,pi,120),r.phase.fit(rec,:),'g')
%%

