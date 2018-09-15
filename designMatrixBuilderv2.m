




function [DmatX, DmatY, motorX] = designMatrixBuilderv2(V,U,designvars,classes,normalization,removal,balance,sampMethod,dropempty)

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
    
    case 'pas'
        % OLD using every single touch;
        %         hx = [V.var.hit{3}' V.var.hit{4}' V.var.hit{5}'];
        %         hy = ones(size(hx,1),1);
        %
        %         FAx = [V.var.FA{3}' V.var.FA{4}' V.var.FA{5}'];
        %         FAy = ones(size(FAx,1),1);
        %
        %         CRx = [V.var.CR{3}' V.var.CR{4}' V.var.CR{5}'];
        %         CRy1 = ones(size(CRx,1),1);
        %
        %         %Repeating matrix for number of touches
        %         hxMotortmp = [];
        %         for i = 1:length(hxMotor)
        %             hxMotortmp = [hxMotortmp ; repmat(hxMotor(i),V.touchNum.hit(i),1)];
        %         end
        %         FAxMotortmp = [];
        %         for i = 1:length(FAxMotor)
        %             FAxMotortmp = [FAxMotortmp ; repmat(FAxMotor(i),V.touchNum.FA(i),1)];
        %         end
        %         CRxMotortmp = [];
        %         for i = 1:length(CRxMotor)
        %             CRxMotortmp = [CRxMotortmp ; repmat(CRxMotor(i),V.touchNum.CR(i),1)];
        %         end
        %
        %         hx = [hx hxMotortmp];
        %         hy(hx(:,3)>=0,:)=[];hx(hx(:,3)>=0,:)=[];
        %         FAx = [FAx FAxMotortmp];
        %         FAy(FAx(:,3)>=0,:)=[];FAx(FAx(:,3)>=0,:)=[];
        %         CRx = [CRx CRxMotortmp] ;
        %         CRy1(CRx(:,3)>=0,:)=[];CRx(CRx(:,3)>=0,:)=[];
        %
        
        %NEW using mean of touches
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxa' hxm' hxp' hxMotor];
        FAx = [FAxa' FAxm' FAxp' FAxMotor];
        CRx = [CRxa' CRxm' CRxp' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
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
        
    case 'touch'
        hx = [V.touchNum.hit']>=1;
        hy = ones(size(hx,1),1);
        FAx = [V.touchNum.FA']>=1;
        FAy = ones(size(FAx,1),1);
        CRx = [V.touchNum.CR']>=1;
        CRy1 = ones(size(CRx,1),1);
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'motor'
        
        hx = [hxMotor hxMotor];
        FAx = [FAxMotor FAxMotor];
        CRx = [CRxMotor CRxMotor] ;
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
    case 'timing'
        hx = [V.var.hit{10}(:,1) hxMotor];
        FAx = [V.var.FA{10}(:,1) FAxMotor];
        CRx = [V.var.CR{10}(:,1) CRxMotor];
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        %           hy = ones(size(hx,1),1);
        %         FAx = [V.var.FA{7}(:,1) V.var.FA{7}(:,2) V.var.FA{7}(:,3) motorPos(FAxmotoridx)'];
        %          FAy = ones(size(FAx,1),1);
        %         CRx = [V.var.CR{7}(:,1) V.var.CR{7}(:,2) V.var.CR{7}(:,3) motorPos(CRxmotoridx)'];
        %          CRy1 = ones(size(CRx,1),1);
        
    case 'decompTime'
                hx = [cellfun(@mean,V.var.hit{7}(:,1:3)) hxMotor];
                FAx = [cellfun(@mean,V.var.FA{7}(:,1:3)) FAxMotor];
                CRx = [cellfun(@mean,V.var.CR{7}(:,1:3)) CRxMotor];
        
        
        %motor build
%         hxnumTouches = cellfun(@numel,V.var.hit{7}(:,1));
%         FAxnumTouches = cellfun(@numel,V.var.FA{7}(:,1));
%         CRxnumTouches = cellfun(@numel,V.var.CR{7}(:,1));
%         
%         hxmotorcomp = [];
%         CRxmotorcomp = [];
%         FAxmotorcomp = [];
%         for b = 1:length(hxnumTouches)
%         hxmotorcomp = [hxmotorcomp ;repmat(hxMotor(b),hxnumTouches(b),1)];
%         end
%         for b = 1:length(FAxnumTouches)
%         FAxmotorcomp = [FAxmotorcomp ;repmat(FAxMotor(b),FAxnumTouches(b),1)];
%         end
%         for b = 1:length(CRxnumTouches)
%         CRxmotorcomp = [CRxmotorcomp ;repmat(CRxMotor(b),CRxnumTouches(b),1)];
%         end
%         
%         
%         hx = [cell2mat(V.var.hit{7}(:,1:3)) hxmotorcomp];
%         FAx = [cell2mat(V.var.FA{7}(:,1:3)) FAxmotorcomp];
%         CRx = [cell2mat(V.var.CR{7}(:,1:3)) CRxmotorcomp];
%         
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
    case 'roll'
        hxmotoridx = V.var.hit{8}(:,2);
        FAxmotoridx = V.var.FA{8}(:,2);
        CRxmotoridx = V.var.CR{8}(:,2);
        
        hx = [V.var.hit{8}(:,1)  motorPos(hxmotoridx)'];
        hy = ones(size(hx,1),1);
        FAx = [V.var.FA{8}(:,1) motorPos(FAxmotoridx)'];
        FAy = ones(size(FAx,1),1);
        CRx = [V.var.CR{8}(:,1) motorPos(CRxmotoridx)'];
        CRy1 = ones(size(CRx,1),1);
        
    case 'ubered'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,1,U,sampMethod);
        %         hx = [hx' V.touchNum.hit' V.licks.oneT.hit' V.licks.oneT.hit'.*V.licks.twoT.hit' V.licks.oneT.hit'.*V.licks.twoT.hit'.*V.licks.threeT.hit'];
        %         mx = [mx' V.touchNum.miss' V.licks.oneT.miss'];
        %         FAx = [FAx' V.touchNum.FA' V.licks.oneT.FA' V.licks.oneT.FA'.*V.licks.twoT.FA' V.licks.oneT.FA'.*V.licks.twoT.FA'.*V.licks.threeT.FA'];
        %         CRx = [CRx' V.touchNum.CR' V.licks.oneT.CR' V.licks.oneT.CR'.*V.licks.twoT.CR' V.licks.oneT.CR'.*V.licks.twoT.CR'.*V.licks.threeT.CR'];
        
        % %          mx = [mx' V.touchNum.miss' V.licks.oneT.miss'];
        
        hx = [hx' V.touchNum.hit' V.var.hit{10}(:,1)];
        FAx = [FAx' V.touchNum.FA' V.var.FA{10}(:,1)];
        CRx = [CRx' V.touchNum.CR' V.var.CR{10}(:,1)];
        
        %          hx = [hx' V.touchNum.hit'];
        %          FAx = [FAx' V.touchNum.FA'];
        %          CRx = [CRx' V.touchNum.CR'];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'theta'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,1,U,sampMethod);
        hx = [hx'];
        FAx = [FAx'];
        CRx = [CRx'];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'maxP'
        hx = [V.var.hit{9}'  hxMotor];
        FAx = [V.var.FA{9}' FAxMotor];
        CRx = [V.var.CR{9}' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
end

switch dropempty
    case 'yes'
        
        
        if strcmp(designvars,'decompTime')
            hx = hx(~isnan(hx(:,1)),:);
            FAx = FAx(~isnan(FAx(:,1)),:);
            CRx = CRx(~isnan(CRx(:,1)),:);
            
        else
            hitcounts = V.touchNum.hit';
            FAcounts = V.touchNum.FA';
            CRcounts = V.touchNum.CR';
            
            hx = hx(hitcounts>0,:);
            FAx = FAx(FAcounts>0,:);
            CRx = CRx(CRcounts>0,:);
            
            droppedTrials = ['hit drops:' num2str(sum(hitcounts==0)) ' FA drops:' num2str(sum(FAcounts==0)) ' CR drops:' num2str(sum(CRcounts==0))]
            
        end
        %Further tossing of NaN values. NaN values will come from
        %touchTimingDecomp that got rid of retraction touches
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        
    case 'no'
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
        
        
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas') ||strcmp(designvars,'decompTime')
            
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