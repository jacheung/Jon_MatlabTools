%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting function for single whisker precision using ROC curves.
% Will bin motor positions into 10 bins and use those for identify FArate
% vs Hitrate. Will only work for continuous uber array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uber_continuousPrecision(U)


popA = cell(1,6);
figure(430);clf
for ms = 1:length(U)
    motors = U{ms}.meta.motorPosition;
    ttype = U{ms}.meta.trialType;
    tcorr = U{ms}.meta.trialCorrect;
    array = [motors' ttype' tcorr'];
    
    [sortedtmp]= binslin(array(:,1),array,'equalE',21,U{ms}.meta.ranges(1),U{ms}.meta.ranges(end));
    [sortedtmp2]= binslin(array(:,1),array,'equalE',11,U{ms}.meta.ranges(1),U{ms}.meta.ranges(end));
    sorted = [sortedtmp2(1:5) ;sortedtmp(10:11); sortedtmp2(6:10)];

    a = [6  5  4 3 2 1];
    b = [7 8 9 10 11 12];
    colors = prism(6);
    for d = 1:length(b)
        onemm = [sorted{a(d)};sorted{b(d)}];
        onemmFA = 1-mean(onemm(onemm(:,2)==0,3));
        onemmHR = mean(onemm(onemm(:,2)==1,3));
        
        popA{d}(ms,:) = [onemmFA onemmHR];
        
%         figure(430);subplot(2,5,ms)
        figure(430);
        hold on; scatter(onemmFA,onemmHR,'markeredgecolor',colors(d,:));
    end
    plot([0 1],[0 1],'-k')
    xlabel('FA rate');ylabel('Hit rate')
%     if ms == 5
%         legend('1mm','2mm','3mm','4mm','5mm')
%     end
end

groupPrecision=cell2mat(cellfun(@mean,popA,'uniformoutput',0)');

    % 95% CI
    for i = 1:size(popA,2)
        x=popA{i};
        SEM = nanstd(x)/sqrt(length(x));               % Standard Error
        ts = tinv([0.05  0.95],length(x)-1);      % T-Score
        cibin(i,:) = abs(ts.*SEM);   %confidence intervals
    end


% figure(431);clf
figure(430);hold on;
for i = 1:length(groupPrecision)
    hold on; 
    errorbar(groupPrecision(i,1),groupPrecision(i,2),cibin(i,2),cibin(i,2),cibin(i,1),cibin(i,1),'x','Color',colors(i,:),'linewidth',1)
    
    scatter(groupPrecision(i,1),groupPrecision(i,2),400,'x','linewidth',2,'markeredgecolor',colors(i,:));
%     h.CData = colors(i,:);
    
    xlabel('FA rate');ylabel('Hit rate')
    legend off
end
set(gca,'ytick',[0 .5 1],'xtick',[0 .5 1],'ylim',[0 1]);
plot([0 1],[0 1],'-k')
axis square
% legend('1mm','2mm','3mm','4mm','5mm','location','southeast')

% set(figure(430), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
% print(figure(430),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\indivPrecision'])
% print(figure(431),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\groupPrecision'])
