




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
        
    case 'kappa'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,6,U,sampMethod);
        
        hx = [hxa' hxMotor];
        FAx = [FAxa' FAxMotor];
        CRx = [CRxa' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
    case 'phase'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxp' hxMotor];
        FAx = [FAxp' FAxMotor];
        CRx = [CRxp' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
    case 'amp'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxa' hxMotor];
        FAx = [FAxa' FAxMotor];
        CRx = [CRxa'  CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
    case 'midpoint'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxm'  hxMotor];
        FAx = [ FAxm'  FAxMotor];
        CRx = [ CRxm'  CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
    case 'mpamp'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        
        hx = [hxm' hxa'  hxMotor];
        FAx = [ FAxm'  FAxa' FAxMotor];
        CRx = [ CRxm' CRxa' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        
    case 'mpphase'
        [hxm,~,FAxm,CRxm] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxm' hxp'  hxMotor];
        FAx = [ FAxm'  FAxp' FAxMotor];
        CRx = [ CRxm' CRxp' CRxMotor];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
    case 'ampphase'
        [hxa,~,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hxa' hxp'  hxMotor];
        FAx = [ FAxa'  FAxp' FAxMotor];
        CRx = [ CRxa' CRxp' CRxMotor];
        
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
        
    case 'countsBinary'
        
        hx = [V.touchNum.hit'];
        hx = hx>0;
        hy = ones(size(hx,1),1);
        FAx = [V.touchNum.FA'];
        FAx = FAx>0;
        FAy = ones(size(FAx,1),1);
        CRx = [V.touchNum.CR'];
        CRx = CRx>0;
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
        hx = [cellfun(@nanmean,V.var.hit{7}(:,1:3)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,1:3)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,1:3)) CRxMotor];
        % DECOMP TIME INDIVIDUALS: :1) timeTotouch from trough :2) theta at trough
        % :3) velocity measured in theta/ms
        
    case 'whiskOnset'
        hx = [cellfun(@nanmean,V.var.hit{7}(:,1)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,1)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,1)) CRxMotor];
        
    case 'timeTotouch'
        hx = [cellfun(@nanmean,V.var.hit{7}(:,1)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,1)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,1)) CRxMotor];
    case 'onsetTheta'
        hx = [cellfun(@nanmean,V.var.hit{7}(:,2)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,2)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,2)) CRxMotor];
    case 'velocity'
        hx = [cellfun(@nanmean,V.var.hit{7}(:,3)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,3)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,3)) CRxMotor];
    case 'Ivelocity'
        
        %         [hxiv,~,FAxiv,CRxiv] = meanVarfinder (V,2,U,sampMethod);
        %         hx = [hxiv' hxMotor];
        %         FAx = [FAxiv' FAxMotor];
        %         CRx = [CRxiv' CRxMotor];
        %
        hx = [cellfun(@nanmean,V.var.hit{7}(:,4)) hxMotor];
        FAx = [cellfun(@nanmean,V.var.FA{7}(:,4)) FAxMotor];
        CRx = [cellfun(@nanmean,V.var.CR{7}(:,4)) CRxMotor];
        
    case 'Ivelocity2'
        
        [hxiv,~,FAxiv,CRxiv] = meanVarfinder (V,2,U,sampMethod);
        hx = [hxiv' hxMotor];
        FAx = [FAxiv' FAxMotor];
        CRx = [CRxiv' CRxMotor];
        %
        
        
        
        
    case 'radialD'
        hx = [V.var.hit{11}(:,1)./33 hxMotor];
        FAx = [V.var.FA{11}(:,1)./33 FAxMotor];
        CRx = [V.var.CR{11}(:,1)./33 CRxMotor];
        
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
        
        [hx,~,FAx,CRx] = meanVarfinder (V,1,U,sampMethod); %theta
        [hxk,~,FAxk,CRxk] = meanVarfinder (V,6,U,sampMethod); %kappa
        %         [hxiv,~,FAxiv,CRxiv] = meanVarfinder (V,2,U,sampMethod);  %vel
        
        hxw = [cellfun(@nanmean,V.var.hit{7}(:,1)) ]; %whisk onset latency
        FAxw = [cellfun(@nanmean,V.var.FA{7}(:,1)) ];
        CRxw = [cellfun(@nanmean,V.var.CR{7}(:,1)) ];
        
        hx = [hxk' V.var.hit{10}(:,1) hxw   V.touchNum.hit' hx']; %10 is cue to touch
        FAx = [FAxk' V.var.FA{10}(:,1) FAxw V.touchNum.FA'  FAx'   ];
        CRx = [CRxk'  V.var.CR{10}(:,1) CRxw  V.touchNum.CR' CRx'];
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    
        
        
    case 'countsmidpoint'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,4,U,sampMethod);
        
        hx = [hx' V.touchNum.hit'];
        FAx = [FAx' V.touchNum.FA'];
        CRx = [CRx' V.touchNum.CR'];
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'countsmpamp'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,4,U,sampMethod);
        [hxa,mxa,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        
        hx = [hx' hxa' V.touchNum.hit'];
        FAx = [FAx' FAxa' V.touchNum.FA'];
        CRx = [CRx' CRxa' V.touchNum.CR'];
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'countsmpphase'
        
        [hx,~,FAx,CRx] = meanVarfinder (V,4,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [hx' hxp' V.touchNum.hit'];
        FAx = [FAx' FAxp' V.touchNum.FA'];
        CRx = [CRx' CRxp' V.touchNum.CR'];
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
        
    case 'countsamp'
        
        [hxa,mxa,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        
        hx = [ hxa' V.touchNum.hit'];
        FAx = [FAxa' V.touchNum.FA'];
        CRx = [CRxa' V.touchNum.CR'];
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
        
    case 'countsampphase'
        
        [hxa,mxa,FAxa,CRxa] = meanVarfinder (V,3,U,sampMethod);
        [hxp,~,FAxp,CRxp] = meanVarfinder (V,5,U,sampMethod);
        
        hx = [ hxa' hxp' V.touchNum.hit'];
        FAx = [FAxa' FAxp' V.touchNum.FA'];
        CRx = [CRxa' CRxp' V.touchNum.CR'];
        
        
        hy = ones(size(hx,1),1);
        FAy = ones(size(FAx,1),1);
        CRy1 = ones(size(CRx,1),1);
        
        hx = [hx hxMotor];
        FAx = [FAx FAxMotor];
        CRx = [CRx CRxMotor] ;
    case 'countsphase'
        
        [hx,mx,FAx,CRx] = meanVarfinder (V,5,U,sampMethod);
        
        
        hx = [hx' V.touchNum.hit'];
        FAx = [FAx' V.touchNum.FA'];
        CRx = [CRx' V.touchNum.CR'];
        
        
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
        
        if strcmp(designvars,'decompTime') || strcmp(designvars,'timeTotouch') ||strcmp(designvars,'onsetTheta') || strcmp(designvars,'velocity') || strcmp(designvars,'Ivelocity') || strcmp(designvars,'whiskOnset') || strcmp(designvars,'ubered')
            
%             hx = hx(~isnan(hx(:,1)),:);
%             FAx = FAx(~isnan(FAx(:,1)),:);
%             CRx = CRx(~isnan(CRx(:,1)),:);
%             
            hx = hx(~sum(isnan(hx),2),:);
            FAx = FAx(~sum(isnan(FAx),2),:);
            CRx = CRx(~sum(isnan(CRx),2),:);
            
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
        
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas') ||strcmp(designvars,'decompTime')
            
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
        
        
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas') ||strcmp(designvars,'decompTime') || strcmp(designvars,'countsMP') || strcmp(designvars,'countsMPamp')
            
            switch normalization
                case 'meanNorm'
                   DmatX = DmatX-mean(DmatX); %mean normalization
                case 'whiten'    
                    DmatX = filex_whiten(DmatX);
            end
            
            
        end
        
    case 'allBehavTypes'
        DmatX=[hx;FAx;CRx]; DmatY = [hy;FAy*2;CRy1*3];
end