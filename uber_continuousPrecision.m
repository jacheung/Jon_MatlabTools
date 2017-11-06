%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting function for single whisker precision using ROC curves.
% Will bin motor positions into 10 bins and use those for identify FArate
% vs Hitrate. Will only work for continuous uber array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uber_continuousPrecision(U)


popA = cell(1,5);
figure(430);clf
for ms = 1:length(U)
    motors = U{ms}.meta.motorPosition;
    ttype = U{ms}.meta.trialType;
    tcorr = U{ms}.meta.trialCorrect;
    array = [motors' ttype' tcorr'];
    
    [sorted]= binslin(array(:,1),array,'equalE',11,U{ms}.meta.ranges(1),U{ms}.meta.ranges(end));
    
    a = [5 4 3 2 1];
    b = [6 7 8 9 10];
    colors = parula(5);
    for d = 1:length(b)
        onemm = [sorted{a(d)};sorted{b(d)}];
        onemmFA = 1-mean(onemm(onemm(:,2)==0,3));
        onemmHR = mean(onemm(onemm(:,2)==1,3));
        
        popA{d}(ms,:) = [onemmFA onemmHR];
        
%         figure(430);subplot(2,5,ms)
        figure(430);
        hold on; h = scatter(onemmFA,onemmHR,'filled');
        h.CData = colors(d,:); alpha(.5)
    end
    plot([0 1],[0 1],'-k')
    xlabel('FA rate');ylabel('Hit rate')
    if ms == 5
        legend('1mm','2mm','3mm','4mm','5mm')
    end
end

groupPrecision=cell2mat(cellfun(@mean,popA,'uniformoutput',0)');
% figure(431);clf
figure(430);hold on;
for i = 1:length(groupPrecision)
    hold on; h=scatter(groupPrecision(i,1),groupPrecision(i,2),400,'x','linewidth',5);
    h.CData = colors(i,:);
    
    xlabel('FA rate');ylabel('Hit rate')
    
end
set(gca,'ytick',[0 .5 1],'xtick',[0 .5 1]);
plot([0 1],[0 1],'-k')
legend('1mm','2mm','3mm','4mm','5mm','location','southeast')

% set(figure(430), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
% print(figure(430),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\indivPrecision'])
% print(figure(431),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\groupPrecision'])
