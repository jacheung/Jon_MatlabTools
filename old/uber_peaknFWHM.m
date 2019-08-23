% window = [8 25];   
basewindow = [-25 0];
figure(570);clf
figure(571);clf
figure(572);clf
recs = [1 3 4 5 9 10 11 12 14];
integWndows=[9 30; 5 30 ; 5 30; 5 30 ; 10 25; 9 20; 10 40; 8 13; 8 25];

for k = 1:length(recs)
    window = integWndows(k,:);
    rec = recs(k);

touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
    V.bounds = [{[-20:1:80]} {[-9750:125:9750]} {[-99.5:5:99.5]} {[-99.5:5:99.5]} {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
    
    
    for d = 1
        binrangetmp = V.bounds{d};
        binspace = (binrangetmp(end)-binrangetmp(1))/(numel(binrangetmp)-1);
        binrange = binrangetmp(1)+binspace/2:binspace:binrangetmp(end)-binspace/2;
        bloop = parula(size(binrange,2));
        
        %Build Variables
        var = squeeze(U{rec}.S_ctk(d,:,:));
        varAtTouchraw = var(touchIdx);
        varAtTouchraw = varAtTouchraw(~isnan(varAtTouchraw));
        touchIdx = touchIdx(~isnan(varAtTouchraw));
        spikesAtTouchraw = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(varAtTouchraw),1)),2);
        spikesPRETouchraw = sum(spikes(repmat(touchIdx,1,numel(basewindow(1):basewindow(2)))+repmat([basewindow(1):basewindow(2)],numel(varAtTouchraw),1)),2);
        
        
        normdSpks = (spikesAtTouchraw-spikesPRETouchraw)./(spikesAtTouchraw+spikesPRETouchraw);
        normdSpks(isnan(normdSpks))=0; %not sure if this correct. setting modidx for 0/0 to 0;
        
   [sorted, ~ ,~]=binslin(varAtTouchraw,spikesAtTouchraw,'equalE',numel(V.bounds{d}),V.bounds{d}(1),V.bounds{d}(end));
     [sortedNormd, ~ ,~]=binslin(varAtTouchraw,normdSpks,'equalE',numel(V.bounds{d}),V.bounds{d}(1),V.bounds{d}(end));
        
%using raw counts 
   means = cellfun(@mean,sorted);
   binrange = binrange(~isnan(means));
   
    smoothedmeans = smooth(means,10);
    smoothedmeans = smoothedmeans(~isnan(means)) ;
    [~,fridx] = max(smoothedmeans);
    [~,fwhmidx] = min(abs(smoothedmeans-(max(smoothedmeans)/2)));
   peak = binrange(fridx);
   fwhm = binrange(fwhmidx);

   
   %using normalized firing rate to baseline fr
   normdmeans = cellfun(@mean,sortedNormd);
   normdmeans = normdmeans(~isnan(normdmeans));
    smoothednormmeans = smooth(normdmeans,10);
    smoothednormmeans = smoothednormmeans(~isnan(normdmeans)) ;
    [~,fridx] = max(smoothednormmeans);
    [~,fwhmidx] = min(abs(smoothednormmeans-(max(smoothednormmeans)/2)));
   normpeak = binrange(fridx);
   normfwhm = binrange(fwhmidx);

   
   % plotting
   figure(570);subplot(5,4,rec)
   plot(binrange,smoothedmeans)
   set(gca,'xlim',[binrange(1) binrange(end)+5],'ylim',[0 max(means)*1.5])
   title([num2str(window(1)) ':' num2str(window(2)) 'ms posttouch'])
   hold on; plot([peak peak],[0 max(means)*1.5],'-.k')
   hold on; plot([fwhm fwhm],[0 max(means)*1.5],'-.r')
    
   
   figure(572);
   hold on; errorbar(peak,k,(peak-fwhm),'horizontal','r')
   hold on; scatter(peak,k,'filled','linewidth',10,'markerfacecolor','black')
   set(gca,'ylim',[0 length(recs)+1],'ytick',[],'yticklabel',[])
   xlabel('Whisker Angle');
   title('Peak Tuning +/- FWHM')
   
   
   figure(571);subplot(5,4,rec)
   plot(binrange,smooth(normdmeans,10),'r')
   set(gca,'xlim',[binrange(1) binrange(end)],'ylim',[0 max(normdmeans)*1.5])
    hold on; plot([normpeak normpeak],[0 max(normdmeans)*1.5],'-.k')
   hold on; plot([normfwhm normfwhm],[0 max(normdmeans)*1.5],'-.r')
    end
end