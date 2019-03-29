 
     designvars = 'ubered';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' 4) 'ubered'
        % 5) 'timing' 6) 'motor' 7) 'decompTime'
        classes = 'lick';
        % 1) 'gonogo' 2) 'FAvsCR' 3) 'lick' 4) allBehavTypes
        
        sample ='bias';
        % 1) 'bias' (takes 70% from each class for train) 2) 'random' just takes
        % random 70% to train
        
        % Only for 'ubered' or 'pas'
        normalization = 'none';
        % 1) 'whiten' 2) 'none';
        
        % Only for 'ubered'
        removal = 'no';
        % 1) 'yes' to remove 0 touch trials and auto classify as CR
        
        balance = 'off';
        
        nanMethod = 'random';
        % 1) random (resample NaN vars from all touches)
        % 2) peakPro (replace NaN using var from max protraction)
        % 3) resampSameD (replace NaN using vas from touches in same trial type)
        
        dropempty = 'no';
        % 1) 'yes' = drop trials with 0 touches
        % 2) 'no' = keep all trials
for rec = 1:length(V)
     
            [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod,dropempty);
            
            lnl{rec} = DmatY==1;
            pthetas{rec} = DmatX(:,1);
            pcounts{rec} = DmatX(:,2);
            pmotorX{rec} = normalize_var(motorX,-1,1)*-1;
end

m=cell2mat(pmotorX');
c = cell2mat(pcounts');
t=cell2mat(pthetas');
l = double(cell2mat(lnl'));

[srt,srtB,bins] =binslin(m,t,'equalE',11,-1,1);
tm = cellfun(@mean,srt);
tms = cellfun(@std,srt);

[srtc,srtBc,binsc] =binslin(m,c,'equalE',11,-1,1);
cm = cellfun(@mean,srtc);

xvals = linspace(-1,1,10);

figure(8);clf;subplot(2,1,1)
errorbar(xvals(1:5),tm(1:5),tms(1:5),'bo-')
hold on;errorbar(xvals(6:10),tm(6:10),tms(6:10),'ro-')
set(gca,'xlim',[-1.1 1.1],'xtick',-1:1:1)

figure(8);subplot(2,1,2)
bar(xvals(1:5),cm(1:5),'b')
hold on;bar(xvals(6:10),cm(6:10),'r')
set(gca,'xlim',[-1.1 1.1],'xtick',-1:1:1,'ytick',0:5:10,'ylim',[0 10])




[srtl,srtbl,binsc] =binslin(c,l,'equalE',16,-.1,15.1);

