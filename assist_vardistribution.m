function [govarx,nogovarx,hit,miss,FA,CR,pregovar,prenogovar] = assist_vardistribution(array,variable,prelixGo,prelixNoGo,varargin);
%Must input array(U{rec}), variable can be 1,2,3,4,5 under varNames, and
%predecision touches. varargin can be a vector with time ranges for
%averaging for the variable (will only work for vel or dkappa). 

%This function will output curves with variable distribution for all
%touches and predecision touches. 
%Top plot will be for the variable distribution sorted by go/nogo.
%Bottom plot will be for the variable distribution sorted by
%Hit(blue)/FA(green)/CR(red)/miss(black)


go=find(array.meta.trialType==1);
nogo=find(array.meta.trialType==0);
govar=cell(1,length(go));pregovar=cell(1,length(go));
nogovar=cell(1,length(nogo));prenogovar=cell(1,length(nogo));


switch variable 
    case 1
    varname = 'Theta';
    case 2
    varname = 'Velocity';
    wndow = varargin{1};%window to look before touch for averaging
    sampspace = [-3000:100:3000];
    case 3
    varname = 'Amplitude';
    case 4
    varname = 'Setpoint';
    case 5
    varname = 'Phase';
    sampspace = [-pi:pi/8:pi];
    case 6
    varname = 'D Kappa';
    wndow = varargin{2};%window to look before touch for averaging
    sampspace = [-.05:.0025:.05];
end

if ~exist('sampspace')
    sampspace = [-50:50];
end
%% Finding variable (setpoint,theta,amp,phase) associated with hit,CR,FA,miss trials

for j = 1:length(go)
    gotouchOnIdx=[];
    gotouchOnIdx = [find(array.S_ctk(9,:,go(j))==1); find(array.S_ctk(12,:,go(j))==1)']; %find indices of touch for current go trial
    if exist('wndow') %when doing averages (for vel or kappa) after touch, need to elim those touches that have indices past 4000 frames
    gotouchOnIdx=gotouchOnIdx(gotouchOnIdx+max(wndow)<3999);
    end
    if variable == 2 || variable ==6
        
        gotmp=zeros(size(wndow,2),length(gotouchOnIdx)); %create matrix of m(time points before/after touch) x n(number of touches)
        pregotmp=zeros(size(wndow,2),length(prelixGo{j}));
        
        for f=1:length(wndow) %fill the above matrix with velocity values around window
            gotmp(f,:)=array.S_ctk(variable,gotouchOnIdx+wndow(f),go(j));
            pregotmp(f,:)=array.S_ctk(variable,prelixGo{j}+wndow(f),go(j));
        end
        
        govar{j}=nanmean(gotmp,1); %get average of all velocity values before touch
        pregovar{j}=nanmean(pregotmp,1);
        
    else
        
        govar{j}=array.S_ctk(variable,gotouchOnIdx,go(j));
        pregovar{j}=array.S_ctk(variable,prelixGo{j},go(j));
    end
end

for k = 1:length(nogo)
    nogotouchOnIdx=[];
    nogotouchOnIdx = [find(array.S_ctk(9,:,nogo(k))==1); find(array.S_ctk(12,:,nogo(k))==1)'];
    if exist('wndow')
    nogotouchOnIdx= nogotouchOnIdx(nogotouchOnIdx+max(wndow)<3999);
    end
    if variable == 2 || variable ==6
        
        nogotmp=zeros(size(wndow,2),length(nogotouchOnIdx));
        prenogotmp=zeros(size(wndow,2),length(prelixNoGo{k}));
        
        for f=1:length(wndow)
            nogotmp(f,:)=array.S_ctk(variable,nogotouchOnIdx+wndow(f),nogo(k));
            prenogotmp(f,:)=array.S_ctk(variable,prelixNoGo{k}+wndow(f),nogo(k));
        end
        
        nogovar{k}=nanmean(nogotmp,1); %get average of velocity before touch
        prenogovar{k} = nanmean(prenogotmp,1);
        
    else
        
        nogovar{k}=array.S_ctk(variable,nogotouchOnIdx,nogo(k));
        prenogovar{k}=array.S_ctk(variable,prelixNoGo{k},nogo(k));
    end
end

%Organizing sorted variables above based on trial correct or not 
hit=pregovar(array.meta.trialCorrect(go)==1);
miss=pregovar(array.meta.trialCorrect(go)==0);
CR=prenogovar(array.meta.trialCorrect(nogo)==1);
FA=prenogovar(array.meta.trialCorrect(nogo)==0);


govarx=cell2mat(govar);nogovarx=cell2mat(nogovar);
pregovarx=cell2mat(pregovar);prenogovarx=cell2mat(prenogovar);
%% Plotting the above but with plot instead of BAR
% TOP plot = go vs nogo (blue vs red) and predecision go vs nogo( cyan vs
% magenta) 
toplot = varargin{3};
switch toplot
    case 'on'
figure(70+variable);clf;subplot(2,1,1);
% plot(histc(govarx,sampspace)/sum(histc(govarx,sampspace)));
hold on;plot(histc(pregovarx,sampspace)/sum(histc(pregovarx,sampspace)),'b');
% hold on;plot(histc(nogovarx,sampspace)/sum(histc(nogovarx,sampspace)),'r');
hold on;plot(histc(prenogovarx,sampspace)/sum(histc(prenogovarx,sampspace)),'r');
switch variable
    case 2 
    set(gca,'xtick',linspace(0,70,7),'xticklabel',[-3000:1000:3000])
    title(['Predecision Go/NoGo Velocity ' num2str(wndow(1)) ' to ' num2str(wndow(end)) 'ms PostTouch Distribution'])
    ylabel('Proportion of Touches')
    case 5 
    set(gca,'xtick',[1:4:17],'xticklabel',{'-pi','-pi/2',0,'pi/2','pi'}) 
    title(['Predecision Go/NoGo ' varname ' Distribution'])
    ylabel('Proportion of Touches')
    case 6 
    set(gca,'xtick',linspace(0,45,5),'xticklabel',[-.05:.025:.05])
    title(['Predecision Go/NoGo Kappa Change ' num2str(wndow(1)) ' to ' num2str(wndow(end)) 'ms PostTouch Distribution'])
    ylabel('Proportion of Touches')
    otherwise
    set(gca,'xtick',[0:25:100],'xticklabel',[-50:25:50],'xlim',[0 100]);
    title(['Predecision Go/NoGo' varname ' Distribution'])
    ylabel('Proportion of Touches')
end


% BOTTOM plot = predecision variables divided between hit/cr/fa/miss 
subplot(2,1,2); %subplot(1,2,2)
plot(histc(cell2mat(hit),sampspace)/sum(histc(cell2mat(hit),sampspace)),'b');
hold on;plot(histc(cell2mat(FA),sampspace)/sum(histc(cell2mat(FA),sampspace)),'g');
hold on;plot(histc(cell2mat(CR),sampspace)/sum(histc(cell2mat(CR),sampspace)),'r');
hold on;plot(histc(cell2mat(miss),sampspace)/sum(histc(cell2mat(miss),sampspace)),'k');
switch variable
    case 2 
    set(gca,'xtick',linspace(0,70,7),'xticklabel',[-3000:1000:3000])
    title(['Go/NoGo Velocity ' num2str(wndow(1)) ' to ' num2str(wndow(end)) 'ms Predecision Touch'])
    ylabel('Proportion of Touches')
    case 5 
    set(gca,'xtick',[1:4:17],'xticklabel',{'-pi','-pi/2',0,'pi/2','pi'}) 
    title(['Predecision ' varname ' @ Touch x Trial Correct'])
    ylabel('Proportion of Touches')
    case 6 
    set(gca,'xtick',linspace(0,45,5),'xticklabel',[-.05:.025:.05])
    title(['Predecision Kappa Change ' num2str(wndow(1)) ' to ' num2str(wndow(end)) 'ms PostTouch Distribution'])
    ylabel('Proportion of Touches')
    otherwise
    set(gca,'xtick',[0:25:100],'xticklabel',[-50:25:50],'xlim',[0 100]);
    title(['Predecision ' varname ' @ Touch x Trial Correct'])
    ylabel('Proportion of Touches')
end
    
end


