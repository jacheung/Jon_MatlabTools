%% find all indices where touch occurs and create window around it. 
touchOnIdx = [find(U{rec}.S_ctk(9,:,:)==1); find(U{rec}.S_ctk(12,:,:)==1)];
touchOffIdx = [find(U{rec}.S_ctk(10,:,:)==1); find(U{rec}.S_ctk(13,:,:)==1)];
spikes = squeeze(U{rec}.R_ntk);
touchOnIdx = touchOnIdx(touchOnIdx<(numel(spikes)-50));
touchOffIdx = touchOffIdx(1:length(touchOnIdx));
touchOnIdx = sort(touchOnIdx,1);
touchOffIdx = sort(touchOffIdx,1);
touchrange=horzcat(touchOnIdx,touchOffIdx);

%% find all params (theta,phase,sp,amp) at touch 

atTouchvar = zeros(length(touchOnIdx),3);

for i = 1:length(touchOnIdx)
    atTouchVar(i,1) = U{rec}.S_ctk(3,touchOnIdx(i)); %amplitude
    atTouchVar(i,2) = U{rec}.S_ctk(4,touchOnIdx(i)); %setpoint
    atTouchVar(i,3) = U{rec}.S_ctk(5,touchOnIdx(i)); %phase
    atTouchVar(i,4) = U{rec}.S_ctk(1,touchOnIdx(i)); %theta
end

%relevancy check (running each IV:DV)
[rho,pval]=corr(atTouchVar,'type','pearson');

figure(10);subplot(1,3,1);scatter(atTouchVar(:,1),atTouchVar(:,4)) %amp x theta
ylabel('theta');xlabel('amplitude')
subplot(1,3,2);scatter(atTouchVar(:,2),atTouchVar(:,4),'r')%setpoint x theta
title('Relevancy Check')
xlabel('setpoint')
subplot(1,3,3);scatter(atTouchVar(:,3),atTouchVar(:,4),'g') %phase x theta
xlabel('phase')

% poor dependency of calculating phase but this may be due to the fact that
% the relationship is NOT linear

%multicolinearity check (running each IV:IV)
figure(3);clf;subplot(1,3,1);scatter(atTouchVar(:,1),atTouchVar(:,2)) %amp x setpoint
    xlabel('amplitude');ylabel('setpoint');
    text(20,6,num2str(rho(1,2)),'FontSize',8,'Color','black')
    text(20,4.5,['p=' num2str(pval(1,2))],'FontSize',8,'Color','black')
subplot(1,3,2);scatter(atTouchVar(:,2),atTouchVar(:,3))%setpoint x phase
    xlabel('setpoint'),ylabel('phase');title('Multicolinearity');
    text(40,-3,num2str(rho(2,3)),'FontSize',8,'Color','black')
    text(40,-3.25,['p=' num2str(pval(2,3))],'FontSize',8,'Color','black')
subplot(1,3,3);scatter(atTouchVar(:,1),atTouchVar(:,3)) %amp x phase
    xlabel('amplitude'),ylabel('phase');
    text(20,-3,num2str(rho(1,3)),'FontSize',8,'Color','black')
    text(20,-3.25,['p=' num2str(pval(1,3))],'FontSize',8,'Color','black')

figure(6);scatter3(atTouchVar(:,1),atTouchVar(:,2),atTouchVar(:,3))
%% OPTICS clustering 
eps=3;
minpts=5;

[m,n] = size(atTouchVar);


clusterd=pdist2(atTouchVar(:,1:2),atTouchVar(:,1:2)); %finding distance b/t all points in amp vs setpoint
%each column = 1 point, each row corresponds to how close that is to another point

CD = zeros(m,1);
RD = ones(m,1)*10*10;
% calculate CORE distance first
for i=1:length(clusterd);
    curr=clusterd(:,i);
    allpts=find(curr<eps);
    if numel(allpts)>=minpts+1%if there are at least min pts
        tmp=sort(curr(allpts));
        CD(i)=tmp(minpts+1);%get CORE distance for that point
        RD(i)=tmp(end); %get reachability D for that point
    else
        CD(i)=100;
        RD(i)=100;
    end
end

% calculate REACHABILITY distance 
order=[];
seeds=[1:m];
ind = 1; %where to start seed

% NEED TO FIGURE OUT HOW TO ORDER

while ~isempty(seeds)
    ob=seeds(ind); %current point you're on       
    seeds(ind)=[]; %eliminate seed you're calculating right now
    order=[order ob]; %add point to seed list
        curr = clusterd(:,ob);
    	mm= find(clusterd(:,ob)<=eps);
        RD(ob)=max(curr(mm)); %finds reachability index for current seed
        
    tmp=find(clusterd(:,order)<=CD(ob));
    ob = tmp(2);
end

%% plotting normalized # of samples on  yaxis and spks/touch on x axis

