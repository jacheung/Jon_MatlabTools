
array=U{1};
follTdist=zeros(1,array.k);
stdfollTdist=zeros(1,array.k);

for i=1:array.k
    touchOnIdx=[];
touchOnIdx = horzcat(find(array.S_ctk(9,:,i)==1),find(array.S_ctk(12,:,i)==1));
foll=[array.whisker.follicleX{i};array.whisker.follicleY{i}];

touchfoll=foll(:,touchOnIdx);
barloc=array.whisker.barPos{i}(2000,2:3);

afollTdist=sqrt(sum((touchfoll-repmat(barloc',1,size(touchfoll,2))).^2));

follTdist(i)=nanmean(afollTdist);
stdfollTdist(i)=nanstd(afollTdist);
end

figure(2);clf;errorbar(array.meta.motorPosition,follTdist/33,stdfollTdist/33,'bo')
xlabel('Motor Position');ylabel('D Follicle-Bar @ Touch')