
for k = 25
    curr = U{k};
    masks = assist_touchmasks(curr);
    
    [~,motorIdx] = sort(curr.meta.motorPosition);
    
    spksraw = squeeze(curr.R_ntk(:,:,:)); 
    amps = squeeze(curr.S_ctk(3,:,:)).*masks.touch; 

    nwper = masks.touch.*(amps<5);
    wper = masks.touch.*(amps>5);
    

    bwmap = [1 1 1;0 0 0];
    cwmap = [0 1 1;1 1 1];
    gwmap = [1 1 1; 0 1 0];
    
    figure(499);clf
    imagesc(masks.touch(:,motorIdx)')
    set(gca,'xtick',[],'ylim',[1 curr.k],'ytick',[])
    colormap(cwmap)
    
    figure(500);clf
    frmotorIdx = flipud(motorIdx');
    for d = 1:curr.k
        st = find(spksraw(:,frmotorIdx(d)));
        if ~isempty(st) 
            hold on;scatter(st,ones(length(st),1)*d,'.k')
        end
    end
    set(gca,'ylim',[0 curr.k+1],'ytick',[])
    
    figure(501);clf
   imagesc(nwper(:,motorIdx)')
    set(gca,'xtick',[],'ylim',[1 curr.k],'ytick',[])
    colormap(bwmap)
    
   figure(502);clf
   imagesc(wper(:,motorIdx)')
    set(gca,'xtick',[],'ylim',[1 curr.k],'ytick',[])
    colormap(gwmap)
    
end