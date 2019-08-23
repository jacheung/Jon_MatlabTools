%need to build A first from varatTouchWindow

touchNames = {'theta','amp','setpoint','phase','pvel','kappa'};
varmod=cell(1,length(touchNames));
[varmod{:}]=deal(zeros(1,length(U)));
win = [6 20];
basewin = [-20 5];
for k = 1:length(touchNames)
    for j = 1:length(U)
        varmean=(mean(A.(touchNames{k}){j}(:,26+win(1):26+win(2)),2)); %parentheses around touchNames{k} = dynamic field names
        %varbase=(mean(A.(touchNames{k}){j}(:,26+basewin(1):26+basewin(2)),2));
        %varbase=mean(varbase);
        %vardiff=varmean-varbase;
        vardiff=varmean;
        binthresh=sum(A.(touchNames{k}){j}(:,78))/numel(A.(touchNames{k}){j}(:,78))*.25;%min number of touches required in each bin
        vardiff=vardiff(A.(touchNames{k}){j}(:,78)>binthresh); %filter out bins with touches less than binthresh
        [~,indmax]=max(abs(vardiff));[~,indmin]=min(abs(vardiff)); %find index of max/min deviation from zero
        varmod{k}(j)=(vardiff(indmax)-vardiff(indmin))/(vardiff(indmax)+vardiff(indmin));
    end
end


tmp=vertcat(varmod{1},varmod{2},varmod{3},varmod{4},varmod{5},varmod{6});
figure(5);imagesc(sortrows(tmp')');colormap(parula);hCBar=colorbar;
set(gca,'yticklabel',touchNames)
title('Modulation Index Sorted by Theta');xlabel('Sorted Cells')