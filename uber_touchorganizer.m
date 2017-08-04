% build object that provides you index of touch, preD,postD
for rec=1:length(U)
DT(rec).touchIdx= [find(U{rec}.S_ctk(9,:,:)==1) ; find(U{rec}.S_ctk(13,:,:)==1)];
lickidx=find(U{rec}.S_ctk(16,:,:)==1)/4000;
firstlix={};
for f = 1:U{rec}.k
    firstlix{f}=find(ceil(lickidx)==f,1);
end
DT(rec).preDtouchIdx
end