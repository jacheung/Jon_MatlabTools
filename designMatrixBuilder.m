



function [DmatX, DmatY, motorX] = designMatrixBuilder(V,U,designvars,classes,normalization,removal,balance,sampMethod)

if strcmp(U.meta.layer,'D')
    U.meta.motorPosition = U.meta.trialType;
    U.meta.motorPosition(U.meta.motorPosition == 0) = -1;
end
motorPos = U.meta.motorPosition;
hxMotor = motorPos(logical(V.trialNums.matrix(1,:)))';
FAxMotor = motorPos(logical(V.trialNums.matrix(3,:)))';
CRxMotor = motorPos(logical(V.trialNums.matrix(4,:)))';
%
switch designvars
    case 'theta'
        hx = [V.var.hit{1}'];
        FAx = [V.var.FA{1}'];
        CRx = [V.var.CR{1}'];
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        
        htmp =[V.touchNum.hit' hxMotor];
        FAtmp = [V.touchNum.FA' FAxMotor];
        CRtmp = [V.touchNum.CR' CRxMotor];
        
        %             htmp =[V.touchNum.hit' ];
        %             FAtmp = [V.touchNum.FA' ];
        %             CRtmp = [V.touchNum.CR' ];
        %
        newMvals = {htmp,FAtmp,CRtmp};
        for d = 1:length(newMvals)%for each trial type
            newmotor = [] ;
            for k=1:length(newMvals{d}) %for trials...
                j = newMvals{d}(k); %for the number of touches within trial
                if ~j==0
                    tmp = repmat(newMvals{d}(k,2),j,1);
                    newmotor = [newmotor; tmp];
                end
            end
            if d == 1
                hx = [hx newmotor];
            elseif d == 2
                FAx = [FAx newmotor];
            elseif d == 3
                CRx = [CRx newmotor];
            end
        end
        
        
    case 'pas'
        %         ntmp=find(V.var.hit{5}<=0);ptmp=find(V.var.hit{5}>0);
        %         hx = [V.var.hit{3}(ntmp)' V.var.hit{4}(ntmp)' V.var.hit{5}(ntmp)'];
        hx = [V.var.hit{3}' V.var.hit{4}' V.var.hit{5}'];
        hy = ones(size(hx,1),1);
        %         Fntmp=find(V.var.FA{5}<=0);Fptmp=find(V.var.FA{5}>0);
        %         FAx = [V.var.FA{3}(Fntmp)' V.var.FA{4}(Fntmp)' V.var.FA{5}(Fntmp)'];
        FAx = [V.var.FA{3}' V.var.FA{4}' V.var.FA{5}'];
        
        FAy = ones(size(FAx,1),1)+1;
        %         Cntmp=find(V.var.CR{5}<=0);Cptmp=find(V.var.CR{5}>0);
        %         CRx = [V.var.CR{3}(Cntmp)' V.var.CR{4}(Cntmp)' V.var.CR{5}(Cntmp)'];
        CRx = [V.var.CR{3}' V.var.CR{4}' V.var.CR{5}'];
        
        CRy1 = ones(size(CRx,1),1)+1;
        
        %Repeating matrix for number of touches 
        hxMotortmp = [];
        for i = 1:length(hxMotor)
            hxMotortmp = [hxMotortmp ; repmat(hxMotor(i),V.touchNum.hit(i),1)];
        end
        FAxMotortmp = [];
        for i = 1:length(FAxMotor)
            FAxMotortmp = [FAxMotortmp ; repmat(FAxMotor(i),V.touchNum.FA(i),1)];
        end 
        CRxMotortmp = [];
        for i = 1:length(CRxMotor)
            CRxMotortmp = [CRxMotortmp ; repmat(CRxMotor(i),V.touchNum.CR(i),1)];
        end
        
        
        hx = [hx hxMotortmp];
        hy(hx(:,3)>=0,:)=[];hx(hx(:,3)>=0,:)=[];
        FAx = [FAx FAxMotortmp];
        FAy(FAx(:,3)>=0,:)=[];FAx(FAx(:,3)>=0,:)=[];
        CRx = [CRx CRxMotortmp] ;
        CRy1(CRx(:,3)>=0,:)=[];CRx(CRx(:,3)>=0,:)=[];
        
    case 'counts'
        
        hx = [V.touchNum.hit'];
        hy = ones(size(hx,1),1);
        FAx = [V.touchNum.FA'];
        FAy = ones(size(FAx,1),1);
        CRx = [V.touchNum.CR'];
        CRy1 = ones(size(CRx,1),1);
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
        
    case 'ubered'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,1,U,sampMethod);
        hx = [hx' V.touchNum.hit' V.licks.oneT.hit' V.licks.oneT.hit'.*V.licks.twoT.hit' V.licks.oneT.hit'.*V.licks.twoT.hit'.*V.licks.threeT.hit'];
        %                 mx = [mx' V.touchNum.miss' V.licks.oneT.miss'];
        FAx = [FAx' V.touchNum.FA' V.licks.oneT.FA' V.licks.oneT.FA'.*V.licks.twoT.FA' V.licks.oneT.FA'.*V.licks.twoT.FA'.*V.licks.threeT.FA'];
        CRx = [CRx' V.touchNum.CR' V.licks.oneT.CR' V.licks.oneT.CR'.*V.licks.twoT.CR' V.licks.oneT.CR'.*V.licks.twoT.CR'.*V.licks.threeT.CR'];
%         
%         hx = [hx' V.touchNum.hit'];
%         %                 mx = [mx' V.touchNum.miss' V.licks.oneT.miss'];
%         FAx = [FAx' V.touchNum.FA'];
%         CRx = [CRx' V.touchNum.CR'];
%         
        
        hy = ones(size(hx,1),1);
        %                 my = ones(size(mx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
end


switch classes
    case 'gonogo'
            DmatX = [hx(:,1:size(hx,2)-1);FAx(:,1:size(FAx,2)-1);CRx(:,1:size(CRx,2)-1)]; 
            DmatY = [hy;FAy.*2;CRy1.*2];
            motorX = [hx(:,size(hx,2));FAx(:,size(FAx,2));CRx(:,size(CRx,2))];
        
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas')
            
            switch removal
                case 'yes'
                    notouchidx= find(DmatX(:,2) == 0);
                    DmatX(notouchidx,:) = [];
                    DmatY(notouchidx,:) = [];
            end
            
            switch normalization
                case 'whiten'
                    DmatX = filex_whiten(DmatX);
            end
            
            
        end
    case 'FAvsCR'
        if strcmp(balance,'on')
            [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
            DmatX = [FAx;CRx]; DmatY = [FAy-1;CRy1];
        else
            DmatX = [FAx;CRx]; DmatY = [FAy-1;CRy1];
        end

    case 'lick'
        if strcmp(balance,'on') %out of date... don't use
            [lix,lixy,nolixx,nolixy] = FACRBalance([hx ; FAx],CRx);
            DmatX = [lix(:,1:size(lix,2)-1);nolixx(:,1:size(lix,2)-1)]; 
            DmatY = [lixy-1;nolixy];
            motorX = [lix(:,size(lix,2));nolixx(:,size(lix,2))];
        else
            DmatX = [hx(:,1:size(hx,2)-1);FAx(:,1:size(FAx,2)-1);CRx(:,1:size(CRx,2)-1)];
            DmatY = [hy;FAy;CRy1.*2];
            motorX = [hx(:,size(hx,2));FAx(:,size(FAx,2));CRx(:,size(CRx,2))];
        end
        
        
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas')
            
            switch removal
                case 'yes'
                    notouchidx= find(DmatX(:,2) == 0);
%                     emptyDmatX = DmatX(notouchidx,:);
%                     emptyDmatY = DmatY(notouchidx,:);
                    DmatX(notouchidx,:) = [];
                    DmatY(notouchidx,:) = [];
            end
            
            switch normalization
                case 'whiten'
                    DmatX = filex_whiten(DmatX);
            end
            
            
        end
        
    case 'allBehavTypes'
        DmatX=[hx;FAx;CRx]; DmatY = [hy;FAy*2;CRy1*3];
end