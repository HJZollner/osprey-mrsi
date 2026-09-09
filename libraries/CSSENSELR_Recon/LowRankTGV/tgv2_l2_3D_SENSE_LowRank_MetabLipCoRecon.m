function [U_rrrc,V_tc, S,LipidData_rrrt,MetabData_rrrt,FreqMap,costFunVal,LipidRM_ff]  = tgv2_l2_3D_SENSE_LowRank_MetabLipCoRecon(Data_ckkkt, alpha0, alpha1,maxits,minits,mrsiReconParams,Threshold)
% Primal dual TGV2 algorithm, as described in the TGV paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute variables / operator for the iterative reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_ckkkt=(single(Data_ckkkt));

[dxm,dym,dzm,dxp,dyp,dzp] = defineDiffOperators();


check_it = mrsiReconParams.LRTGVModelParams.check_it;
Plot_it= mrsiReconParams.LRTGVModelParams.Plot_it;
CorrB0Map_it=100;mrsiReconParams.LRTGVModelParams.CorrB0Map_it;
CorrB0Map_Maxcount=5;mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount;
Orthogonalize_it=mrsiReconParams.LRTGVModelParams.Orthogonalize_it;
reduction =mrsiReconParams.LRTGVModelParams.reduction;

min_taup=mrsiReconParams.LRTGVModelParams.min_taup;
max_taup=mrsiReconParams.LRTGVModelParams.max_taup;
DualPrimalTauFact=mrsiReconParams.LRTGVModelParams.DualPrimalTauFact;
Momentum=mrsiReconParams.LRTGVModelParams.Momentum;

Nlip=32; % for plotting only

[~,MaxLipPPM_pt]=min(abs(mrsiReconParams.LipidMaxPPM  - mrsiReconParams.ppm));
[~,MinLipPPM_pt]=min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));
FreqOI=mrsiReconParams.MinPPM_pt:mrsiReconParams.MaxPPM_pt;


NbTShort=numel(FreqOI);

LipExp=mrsiReconParams.LipidRemovalParams.LipidRemExpo;
TolLipRMDamping=mrsiReconParams.LipidRemovalParams.TolLipRMDamping;

fprintf(['Lipid suppression with tolerance : ',num2str(TolLipRMDamping),' and exponent : ',num2str(LipExp),'\n']);


[ NbCoil M N Slc NbT] = size(Data_ckkkt); % numSamplesOnSpoke, numSamplesOnSpoke, nCh
if (Slc==1)
    Is2D=1; %Dataset is 2D
else
    Is2D=0;
end
SBFreqMax=NbT-mrsiReconParams.MaxPPM_pt;
NbComp = mrsiReconParams.modelOrder; % numSamplesOnSpoke, numSamplesOnSpoke, nCh

SizeVol = [M N Slc NbComp];
UIndcs  = repmat({':'}, [1, numel(SizeVol)]);
numSpatialPts=M*N*Slc;
FreqMap=single(mrsiReconParams.WaterFreqMap) ;
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;

ShortTime=[0 :(NbT-1)]/Fs;


%for Forward Transform

Time_111t=single(reshape(([0 :(NbT-1)]'/Fs),[1 1 1 NbT]));

SignalMask_rrr=(single(mrsiReconParams.ImMask));
SignalMask_rrr1=(single(reshape(mrsiReconParams.ImMask,[size(mrsiReconParams.ImMask) 1])));
SignalMask_1rrr1=(single(reshape(mrsiReconParams.ImMask,[1 size(mrsiReconParams.ImMask) 1])));

FreqRange=([0 :(NbT-1)]-round(NbT/2));

BrainMask=mrsiReconParams.BrainMask;

SENSE_crrr1=(single(reshape(mrsiReconParams.SENSE,[size(mrsiReconParams.SENSE) 1])));
SENSE_crrr1=(single(SENSE_crrr1.* reshape(SignalMask_rrr1,[1 size(SignalMask_rrr1)])));

HannF_1kkk1=(single(reshape(mrsiReconParams.HKernel,[1 size(mrsiReconParams.HKernel) 1])));
kmask_1kkk1=(single(reshape(mrsiReconParams.kmask,[1 size(mrsiReconParams.kmask) 1])));

kmask=single(mrsiReconParams.kmask);
ImMask=single(mrsiReconParams.ImMask);

SkMask_rrrt = single(repmat(mrsiReconParams.SkMask,[1 1 1 NbT]));
ImMask_rrrt = single(repmat(mrsiReconParams.ImMask,[1 1 1 NbT]));


StepNormGrad = [];
StepDiff=Threshold;
StepDiffu = [];
StepDiffr = [];
StepDiffxi = [];
StepDiffq = [ ];
StepDiffp = [];
StepDiffww = [ ];
StepDiffdivp = [];

tau_p = max_taup;%1/128;%64;%%1/64;%1/16; %the lowest the most stable against singular point
tau_d= tau_p*DualPrimalTauFact;%1/64;%1/8;


PrevNormGrad=1E16;
NormGradWentDown=0;
NbDivergSteps=0;

RelDiffSq_V=1;
RelDiffSq_U=1;

t_old = 2;
PrevNormGrad=0;
step_noConv=0;

%Temporary variable for Forward / Adjoint transform
MGrad_rrrt = zeros([M,N,Slc,NbT],class(Data_ckkkt));
Grad_rrrt = zeros([M,N,Slc,NbT],class(Data_ckkkt));
MGradMetab_rrrt = zeros([M,N,Slc,NbT],class(Data_ckkkt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the lipid recon without LR and TGV but Freq Map adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("Perform the lipid recon ...\n");
PrevNormGrad=0;tau=1/8;
UpdatePreReconData_rrrt=0*MGrad_rrrt;
RelDiffPreRecon=1;
CorrB0Map_count=0;

k=-1; %iteration index
CorrB0Map_count=0;

%store variable to GPU
Vars=whos;
thisgpu =SelectFreeGPU(1.5*sum([Vars.bytes]));

if thisgpu>0
	fprintf("Running Recon on GPU nb: "+num2str(thisgpu)+" ...\n");
	
	ForwMaskedDFT_rk=ForwDFT3_Masked_Operator([M N Slc],kmask,ImMask);
	CombMaskedDFT_rr=ForwMaskedDFT_rk*diag(mrsiReconParams.HKernel(kmask>0))*ForwMaskedDFT_rk';
	AdjData_ctr=reshape(permute(Data_ckkkt,[1,5,2,3,4]),[NbCoil*NbT M*N*Slc]);
	AdjData_ctr=AdjData_ctr(:,kmask>0)*diag(mrsiReconParams.HKernel(kmask>0))*ForwMaskedDFT_rk';
	
	AdjRecon_ctrrr=zeros([NbCoil*NbT M*N*Slc]);
	AdjRecon_ctrrr(:,ImMask>0) = AdjData_ctr;
	AdjRecon_ctrrr = reshape(AdjRecon_ctrrr,[NbCoil NbT M N Slc]);
	PreReconData_rrrt =sqz(sum(conj(SENSE_crrr1).*permute(AdjRecon_ctrrr,[1,3,4,5,2]),1));
	clear Data_ckkkt AdjRecon_ctrrr
		
	VarNames=who;
	for VN=1:numel(VarNames)
		 eval([VarNames{VN}, ' = transferToGPU(',VarNames{VN},');']);
	end
	
else
	fprintf("Running Recon on CPU ...\n");
	PreReconData_rrrt =   sqz(sum(conj(SENSE_crrr1).*ifft(ifft(ifft(Data_ckkkt,[],2),[],3),[],4),1));
end


while RelDiffPreRecon>2E-3 & k<200;
    k=k+1;

    if thisgpu>0 % Use a DFT operator that is fatser on GPU
    
	    temp_ctr = single(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*SENSE_crrr1.*reshape(PreReconData_rrrt,[1 M N Slc NbT]) ; %crrrt
	    temp_ctr = reshape(permute(temp_ctr,[1,5,2,3,4]),[NbCoil*NbT M*N*Slc]);
	    temp_ctr(:,ImMask>0)=temp_ctr(:,ImMask>0)*CombMaskedDFT_rr - AdjData_ctr;

	    temp_crrrt=permute(reshape(temp_ctr,[NbCoil NbT M N Slc]),[1,3,4,5,2]);
	    MGrad_rrrt = MGrad_rrrt + tau_d*SignalMask_rrr1.*reshape(sum( conj(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*conj(SENSE_crrr1).*temp_crrrt ,1),[M N Slc NbT]);% r-r-r-t
	    
    else % fft is faster on CPU
            temp_crrrt = single(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*SENSE_crrr1.*reshape(PreReconData_rrrt,[1 M N Slc NbT]) ; %crrrt
	    temp_crrrt = (fft(fft(fft(temp_crrrt,[],2),[],3),[],4) - Data_ckkkt).*kmask_1kkk1.*HannF_1kkk1;%c-k-k-k-t
	    MGrad_rrrt = MGrad_rrrt + tau_d*SignalMask_rrr1.*reshape(sum( conj(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*conj(SENSE_crrr1).*ifft(ifft(ifft(temp_crrrt,[],2),[],3),[],4) ,1),[M N Slc NbT]);% r-r-r-t
	    
    end
    
    MGrad_rrrt = MGrad_rrrt/(1+tau_d);

    UpdatePreReconData_rrrt = - tau*(MGrad_rrrt ) + Momentum*UpdatePreReconData_rrrt; %lines 8
    PreReconData_rrrt = PreReconData_rrrt + UpdatePreReconData_rrrt;
    
    RelGradNormLip= norm(MGrad_rrrt(:))/norm(PreReconData_rrrt(:));

    NormGrad=norm(MGrad_rrrt(:));
    if k>2 && NormGrad<PrevNormGrad
        NormGradWentDown=1;
        tau = tau*1.10;
        if tau>0.5; tau=0.5;end
    elseif NormGradWentDown==1
        tau=tau*0.85;
        if tau<(1/32);tau=(1/32);end
    end
    tau_d=tau*DualPrimalTauFact;
    EneregyPreReconData = sum(abs(PreReconData_rrrt(:)).^2);
    RelDiffPreRecon = norm(UpdatePreReconData_rrrt(SkMask_rrrt>0))/norm(PreReconData_rrrt(SkMask_rrrt>0));
    PrevNormGrad=NormGrad;

  %  fprintf(['Lipid Recon, k =',num2str(k),', EneregyPreReconData=',num2str(EneregyPreReconData),', RelDiffPreRecon=',num2str(RelDiffPreRecon),', tau=',num2str(tau),'\n']);
end

LipidData_rrrt = SkMask_rrrt.*PreReconData_rrrt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute the initial condition for Metab Components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Initial Solution
PreReconData_rrrf = fft(PreReconData_rrrt,[],4);

Lipid_stack_rf=fft(reshape(LipidData_rrrt(SkMask_rrrt>0),[],NbT),[],2);

LipReg=1;
diagLRM=1;
LipExp=1;
while diagLRM>TolLipRMDamping
	if LipReg*1.2<5E5
	    LipReg=LipReg*1.2;
    	else
    	    LipExp=LipExp*1.05;
    	end
    LipidRM_ff = (inv( eye(NbT) + LipReg * Lipid_stack_rf'*(Lipid_stack_rf)/norm(Lipid_stack_rf)^2 ))^LipExp;
    diagLRM=mean(abs(diag(LipidRM_ff)));
end



if mrsiReconParams.SideBandRMRank>0
	fprintf(['Performing sideband removal in reconstruction with rank ',num2str(mrsiReconParams.SideBandRMRank),'...\n']);
	DownField_rf = reshape(PreReconData_rrrf,[],NbT);
	DownField_rf = DownField_rf(BrainMask(:)>0,:);
		
	HalfDownField_rt= ifft(DownField_rf(:,round(end/2):end),[],2);
	SideBandsSignal_rf =  [-fft(conj(HalfDownField_rt),[],2), fft((HalfDownField_rt),[],2)];	
		
	SideBandsSignal_rt = ifft(SideBandsSignal_rf,[],2);
	SideBandsSignal_rf = fft(SideBandsSignal_rt(:,1:NbT),[],2);% truncate the time serie to the right length (is 1 too long if NbT is odd)
	
	[~,S_SB,SB_V_fc] = svd(SideBandsSignal_rf,0);
	clear DownField_rf HalfDownField_rt SideBandsSignal_rf SideBandsSignal_rt

	SB_V_fc=SB_V_fc(:,1:mrsiReconParams.SideBandRMRank);
	SBRMOp_ff=eye(NbT)-(SB_V_fc)*SB_V_fc';
else
	fprintf(['No sideband removal in reconstruction with rank ',num2str(mrsiReconParams.SideBandRMRank),'...\n']);
	SBRMOp_ff=eye(NbT);
end

NormGradLWentDown=0;
MetabData_rrrt= ifft(reshape(reshape(PreReconData_rrrf,[],NbT)*LipidRM_ff*SBRMOp_ff,size(PreReconData_rrrt)),[],4);
MetabData_rrrf =fft(MetabData_rrrt,[],4);

clear UpdatePreReconData_rrrt PreReconData_rrrt PreReconData_rrrf;

[Uorig,Sorig,Vorig] = svd(reshape( BrainMask.*ifft(MetabData_rrrf(:,:,:,FreqOI),[],4),[M*N*Slc NbTShort ]),0);
clear MetabData_rrrf;
V_tc=Vorig(:,1:mrsiReconParams.modelOrder);
S=Sorig(1:mrsiReconParams.modelOrder,1:mrsiReconParams.modelOrder);
U_rrrc=Uorig(:,1:mrsiReconParams.modelOrder);
U_rrrc = BrainMask.*reshape(U_rrrc,[M N Slc mrsiReconParams.modelOrder]);

U_rrrc=(single(U_rrrc));
V_tc=(single(V_tc));
S=(single(S));

[Uorig,Sorig,Vorig] = svd(reshape( SkMask_rrrt.*LipidData_rrrt ,[M*N*Slc NbT ]),0);
Vlip_tc=Vorig(:,1:Nlip);
Slip=Sorig(1:Nlip,1:Nlip);
Ulip_rrrc=Uorig(:,1:Nlip);
Ulip_rrrc = SignalMask_rrr1.*reshape(Ulip_rrrc,[M N Slc Nlip]);

VisualizeTGV( U_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_initial_Recon_Metabs']);
VisualizeSpectral( V_tc,S, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_initial_AfterLipPreRecon_Metabs']);
VisualizeTGV( Ulip_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_initial_Recon_Lipids',]);
VisualizeSpectral( Vlip_tc,Slip, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_initial_AfterLipPreRecon_Lipids']);
if mrsiReconParams.SideBandRMRank>0
	VisualizeSpectral( ifft(SB_V_fc,[],1),S_SB, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_initial_SideBands']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Metabolite LR & TGV Recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_rrrc=U_rrrc.*SignalMask_rrr1;
u = U_rrrc;
V_proj    = V_tc;

p = zeros([SizeVol,3],class(MetabData_rrrt));
Du = zeros([SizeVol,3],class(MetabData_rrrt));
q = zeros([SizeVol,6],class(MetabData_rrrt));
Dq = zeros([SizeVol,3],class(MetabData_rrrt));
xi = zeros([SizeVol,3],class(MetabData_rrrt));
Dxi = zeros([SizeVol,6],class(MetabData_rrrt));
xi_ = xi;

% Compute the initial data consistency value
dataConsistencyCost=0;
costFunVal = dataConsistencyCost;
TGVCost0Val = 0;
TGVCost1Val = 0;

alpha00 = alpha0/reduction;
alpha10 = alpha1/reduction;
alpha002 = alpha0*reduction;
alpha102 = alpha1*reduction;
alpha01 = alpha0;
alpha11 = alpha1;

UpdateLipidData_rrrt = LipidData_rrrt;

MetabDataLR_rrrt =0*MetabData_rrrt;
UpdateMetabData_rrrt=0*MetabData_rrrt;

MGrad_rrrt=0*MGrad_rrrt;
RGradLipRM_rrrt=0*MGrad_rrrt;
k=-1; %iteration index
CorrB0Map_count=0;
ExtrLipRMinProgress=0;

%store variable to GPU

if thisgpu>0
	VarNames=who;
	for VN=1:numel(VarNames)
		 eval([VarNames{VN}, ' = transferToGPU(',VarNames{VN},');']);
	end
end

while k<maxits %(abs(StepDiff(end))>Threshold  & (k<maxits) ) | k<minits;

    k=k+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL CONVERGENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % update alpha's

    if k<=(minits/2)   
        
	alpha0 = exp(k/(minits/2)*log(alpha01) + ((minits/2)-k)/(minits/2)*log(alpha002));
        alpha1 = exp(k/(minits/2)*log(alpha11) + ((minits/2)-k)/(minits/2)*log(alpha102));

    else

        alpha0 = alpha01;
        alpha1 = alpha11;
    end

    if k > minits
        k_ind=(k-minits)/(maxits-minits); %from 0 to 1
        max_SpectStep=mrsiReconParams.LRTGVModelParams.max_SpectStep;
        max_taup = mrsiReconParams.LRTGVModelParams.max_taup;
        max_taup = min_taup + (max_taup-min_taup)*exp(-k_ind*10);
        
    end

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DUAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Is2D==0 %3D case

        % gradient
        Du(UIndcs{:},1)  = dxp(U_rrrc);
        Du(UIndcs{:},2)  = dyp(U_rrrc);
        Du(UIndcs{:},3) = dzp(U_rrrc);
       
	p = p - tau_d*(Du + xi_);

        % projection

	for comp=1:NbComp
            absp = sqrt(abs(p(UIndcs{1:(end-1)},comp,1)).^2 + abs(p(UIndcs{1:(end-1)},comp,2)).^2  + abs(p(UIndcs{1:(end-1)},comp,3)).^2);
            denom = max(1,absp/(alpha1*max(diag(S))/S(comp,comp)));
            p(UIndcs{1:(end-1)},comp,:) = p(UIndcs{1:(end-1)},comp,:)./denom;

       end

        % symmetrized gradient
       
        Dxi(UIndcs{:},1) = dxm(xi_(UIndcs{:},1));
         Dxi(UIndcs{:},2)  = dym(xi_(UIndcs{:},2));
         Dxi(UIndcs{:},3)  = dzm(xi_(UIndcs{:},3));
        Dxi(UIndcs{:},4)  = (dxm(xi_(UIndcs{:},2)) + dym(xi_(UIndcs{:},1)))/2;
         Dxi(UIndcs{:},5)  = (dxm(xi_(UIndcs{:},3)) + dzm(xi_(UIndcs{:},1)))/2;
        Dxi(UIndcs{:},6)  = (dym(xi_(UIndcs{:},3)) + dzm(xi_(UIndcs{:},2)))/2;
       
	q = q - tau_d*Dxi;

        % projection
        
          for comp=1:NbComp
            absq = sqrt(abs(q(UIndcs{1:(end-1)},comp,1)).^2 + abs(q(UIndcs{1:(end-1)},comp,2)).^2 + abs(q(UIndcs{1:(end-1)},comp,3)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,4)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,5)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,6)).^2);
            denom = max(1,absq/(alpha0*max(diag(S))/S(comp,comp)));
            q(UIndcs{1:(end-1)},comp,:) = q(UIndcs{1:(end-1)},comp,:)./denom;
        end
        
    else	 %2D case

        % gradient
        ux = dxp(U_rrrc);
        uy = dyp(U_rrrc);

        p(UIndcs{:},1) = p(UIndcs{:},1) - tau_d*(ux + xi_(UIndcs{:},1));
        p(UIndcs{:},2) = p(UIndcs{:},2) - tau_d*(uy + xi_(UIndcs{:},2));

        % projection

        for comp=1:NbComp
            absp = sqrt(abs(p(UIndcs{1:(end-1)},comp,1)).^2 + abs(p(UIndcs{1:(end-1)},comp,2)).^2 );
            denom = max(1,absp/(alpha1*max(diag(S))/S(comp,comp)));
            % denom = max(1,absp/(alpha1));
            p(UIndcs{1:(end-1)},comp,1) = p(UIndcs{1:(end-1)},comp,1)./denom;
            p(UIndcs{1:(end-1)},comp,2) = p(UIndcs{1:(end-1)},comp,2)./denom;
        end


        % symmetrized gradient
        gradxi1 = dxm(xi_(UIndcs{:},1));
        gradxi2 = dym(xi_(UIndcs{:},2));
        gradxi3 = (dym(xi_(UIndcs{:},1)) + dxm(xi_(UIndcs{:},2)))/2;

        q(UIndcs{:},1) = q(UIndcs{:},1) - tau_d*gradxi1; % line
        q(UIndcs{:},2) = q(UIndcs{:},2) - tau_d*gradxi2;
        q(UIndcs{:},3) = q(UIndcs{:},3) - tau_d*gradxi3;

        % projection

        for comp=1:NbComp
            absq = sqrt(abs(q(UIndcs{1:(end-1)},comp,1)).^2 + abs(q(UIndcs{1:(end-1)},comp,2)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,3)).^2);
            denom = max(1,absq/(alpha0*max(diag(S))/S(comp,comp)));
            % denom = max(1,absq/(alpha0));
            q(UIndcs{1:(end-1)},comp,1) = q(UIndcs{1:(end-1)},comp,1)./denom;
            q(UIndcs{1:(end-1)},comp,2) = q(UIndcs{1:(end-1)},comp,2)./denom;
            q(UIndcs{1:(end-1)},comp,3) = q(UIndcs{1:(end-1)},comp,3)./denom;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIMAL SPATIAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
    MetabDataLR_rrrf = fft(MetabData_rrrt,[],4);
    MetabDataLRShort_rrrf= fft(formTensorProduct(U_rrrc, V_tc*S),[],4);
    MetabDataLR_rrrf(:,:,:,FreqOI) = MetabDataLRShort_rrrf;
    		
    if thisgpu>0 % Use a DFT operator that is fatser on GPU
	    temp_ctr = single(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*SENSE_crrr1.*reshape(ifft(MetabDataLR_rrrf,[],4) + LipidData_rrrt ,[1 M N Slc NbT]) ; %crrrt
	    temp_ctr = reshape(permute(temp_ctr,[1,5,2,3,4]),[NbCoil*NbT M*N*Slc]);
	    temp_ctr(:,ImMask>0)=temp_ctr(:,ImMask>0)*CombMaskedDFT_rr-AdjData_ctr;

	    temp_crrrt=permute(reshape(temp_ctr,[NbCoil NbT M N Slc]),[1,3,4,5,2]);
	    Grad_rrrt = SignalMask_rrr1.*reshape(sum( conj(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*conj(SENSE_crrr1).*temp_crrrt ,1),[M N Slc NbT]);% r-r-r-t
    else % fft is faster on CPU	   

	    temp_crrrt = single(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*SENSE_crrr1.*reshape(ifft(MetabDataLR_rrrf,[],4) + LipidData_rrrt ,[1 M N Slc NbT]) ;
	    temp_crrrt = (fft(fft(fft(temp_crrrt,[],2),[],3),[],4) - Data_ckkkt).*kmask_1kkk1.*HannF_1kkk1;%c-k-k-k-t 
	    Grad_rrrt = SignalMask_rrr1.*reshape(sum( conj(reshape(exp(2*pi*1i*Time_111t.*FreqMap),[1 M N Slc NbT])).*conj(SENSE_crrr1).*ifft(ifft(ifft(temp_crrrt,[],2),[],3),[],4) ,1),[M N Slc NbT]);% r-r-r-t
    end
   
    MGrad_rrrt = (MGrad_rrrt+tau_d*Grad_rrrt)/(1+tau_d);
    
 if ( mod(k,5)==1)

    	Lipid_stack_rf = fft(reshape(LipidData_rrrt(SkMask_rrrt>0),[],NbT),[],2);

	LipReg=LipReg/10;
	diagLRM=1;
	LipExp=1;
	while diagLRM>TolLipRMDamping
		if LipReg*1.2<5E5
		    LipReg=LipReg*1.2;
	    	else
	    	    LipExp=LipExp*1.05;
	    	end
	    LipidRM_ff = (inv( eye(NbT) + LipReg * Lipid_stack_rf'*(Lipid_stack_rf)/norm(Lipid_stack_rf)^2 ))^LipExp;
	    diagLRM=mean(abs(diag(LipidRM_ff)));
	end

 end   
 	


     MGradMetab_rrrt =( MGradMetab_rrrt + tau_d*ifft(reshape(reshape(fft(Grad_rrrt,[],4),[],NbT)*LipidRM_ff*SBRMOp_ff,size(Grad_rrrt)),[],4) ) /(1+tau_d);

if (  k>50  && k<(minits) )	   
	    RGradLipRM_rrrt =  RGradLipRM_rrrt + tau_d*ifft( MetabDataLR_rrrf - reshape(reshape(MetabDataLR_rrrf,[],NbT)*LipidRM_ff*SBRMOp_ff,size(MetabDataLR_rrrf)),[],4) ;	   
	    RGradLipRM_rrrt = RGradLipRM_rrrt/(1+tau_d); 
else
	    RGradLipRM_rrrt=0*MetabDataLR_rrrf; 
end

    
    UpdateLipidData_rrrt = UpdateLipidData_rrrt - tau_p*SkMask_rrrt.*(MGrad_rrrt);  
    LipidData_rrrt = UpdateLipidData_rrrt + Momentum*(UpdateLipidData_rrrt - LipidData_rrrt);
    
    UpdateMetabData_rrrt = UpdateMetabData_rrrt - tau_p*SignalMask_rrr1.*(MGradMetab_rrrt + RGradLipRM_rrrt);  
    MetabData_rrrt = UpdateMetabData_rrrt + Momentum*(UpdateMetabData_rrrt - MetabData_rrrt);     
   
    RelDiffMetab = norm(UpdateMetabData_rrrt(ImMask_rrrt>0) )/norm(MetabData_rrrt(ImMask_rrrt>0));
    RelDiffLipid = norm(UpdateLipidData_rrrt(SkMask_rrrt>0))/norm(LipidData_rrrt(SkMask_rrrt>0)); 
            
    RGradMetabShort_rrrt=fft(MGradMetab_rrrt + RGradLipRM_rrrt,[],4);
    RGradMetabShort_rrrt=ifft(RGradMetabShort_rrrt(:,:,:,FreqOI),[],4);     
             
    ww = formTensorProduct(RGradMetabShort_rrrt ,(V_tc*diag(1./diag(S)))' ) ;%r-r-r-c
   	
    if Is2D==0 %3D case
        % divergence
        divp = dxm(p(UIndcs{:},1)) + dym(p(UIndcs{:},2)) + dzm(p(UIndcs{:},3));

        u = SignalMask_rrr1.*(u - tau_p*(ww + divp)); %lines 8

        % divergence
        
        Dq(UIndcs{:},1) = dxp(q(UIndcs{:},1)) + dyp(q(UIndcs{:},4)) + dzp(q(UIndcs{:},5));
        Dq(UIndcs{:},2) = dxp(q(UIndcs{:},4)) + dyp(q(UIndcs{:},2)) + dzp(q(UIndcs{:},6));
        Dq(UIndcs{:},3) = dxp(q(UIndcs{:},5)) + dyp(q(UIndcs{:},6)) + dzp(q(UIndcs{:},3));
             
        xi = xi - tau_p*(Dq - p); %line 11
    else %2D case
        % divergence
        divp = dxm(p(UIndcs{:},1)) + dym(p(UIndcs{:},2));

        u = u - tau_p*(ww + divp); %lines 8
	
        % divergence
        divq1 = dxp(q(UIndcs{:},1)) + dyp(q(UIndcs{:},3));
        divq2 = dxp(q(UIndcs{:},3)) + dyp(q(UIndcs{:},2));

        xi(UIndcs{:},1) = xi(UIndcs{:},1) - tau_p*(divq1 - p(UIndcs{:},1));%line 11
        xi(UIndcs{:},2) = xi(UIndcs{:},2) - tau_p*(divq2 - p(UIndcs{:},2));%line 11

    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIARY SPATIAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if  mod(k+1,Orthogonalize_it)==0
        U_rrrc = reshape(OrthogonalizeComponents(reshape(U_rrrc,[size(U_rrrc,1)*size(U_rrrc,2)*size(U_rrrc,3) , size(U_rrrc,4)])),size(U_rrrc));
        u = reshape(OrthogonalizeComponents(reshape(u,[size(u,1)*size(u,2)*size(u,3) , size(u,4)])),size(u));
    end


    for c=1:NbComp
        NormUrrc=sqrt(sum(sum(sum(abs(u(:,:,:,c)).^2,1),2),3));
        S(c,c)=S(c,c)*NormUrrc;
        U_rrrc(:,:,:,c) = U_rrrc(:,:,:,c)/NormUrrc;
        u(:,:,:,c) = u(:,:,:,c)/NormUrrc;
    end


    U_rrrc = u + Momentum*(u - U_rrrc);
    
    xi_ = xi + Momentum*(xi - xi_);


    RelDiffSq_U = norm(u(:))^2/norm(U_rrrc(:))^2 + norm(xi_(:)-xi(:))^2/norm(xi(:))^2;
    norm_p=norm(p(:),1);
    norm_q=norm(q(:),1);
    RelDiffSq_old=RelDiffSq_U;

	NormGrad=norm(RGradMetabShort_rrrt(:));
    if k>10 && NormGrad<PrevNormGrad
        NormGradWentDown=1;
        tau_p = tau_p*1.10;
        if tau_p>max_taup; tau_p=max_taup;end
    elseif NormGradWentDown==1
        tau_p=tau_p*0.85;%0.2;
        if tau_p<min_taup;tau_p=min_taup;end
    end
    
    tau_d=tau_p*DualPrimalTauFact;
    PrevNormGrad=NormGrad;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIMAL SPECTRAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gradStep = ((reshape(U_rrrc, [], NbComp)*diag(1./diag(S)))' * ( reshape(RGradMetabShort_rrrt , [numSpatialPts, size(RGradMetabShort_rrrt,4)]) ) )' ; %t -c
    if norm(gradStep)>0
        [ gradStep] = OrthogonalizeComponents(gradStep);
    end
    V_proj   = V_proj -tau_p * gradStep ;%

    for c=1:size(V_tc,2)
        V_proj(:,c) = V_proj(:,c) ./ norm(V_proj(:,c));
    end

    if k>0
        RelDiffSq_V=norm(V_proj(:))^2/norm(V_tc(:))^2;
    end
    

    V_tc = V_proj + Momentum*( V_proj - V_tc);
    
    for c=1:size(V_tc,2)
        V_tc(:,c) = V_tc(:,c) ./ norm(V_tc(:,c));
    end

    if  mod(k+1,Orthogonalize_it)==0
        V_tc = OrthogonalizeComponents(V_tc);
        V_proj= OrthogonalizeComponents(V_proj);

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FREQUENCY MAP DYNAMIC CORRECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (round(mod(k,CorrB0Map_it)) == (CorrB0Map_it-1) ) & ( CorrB0Map_count<=CorrB0Map_Maxcount ) & (k<minits)
        CorrB0Map_count=CorrB0Map_count+1;
        FreqRange=mrsiReconParams.WaterRefFreqSearchRange;
        FreqPrec=0.5;

        %Metabolites
        A=8;
        RefTimeSerie=conj(sum((V_tc(:,1:A)*S(1:A,1:A)),2));
        if (CorrB0Map_count==1)
            InitMetabRefTimeSerie=RefTimeSerie;
        end


        Fs = mrsiReconParams.mrProt.samplerate*numel(RefTimeSerie)/mrsiReconParams.mrProt.VSize;
        Time = ([0 :(numel(RefTimeSerie)-1)]'/(Fs));
        [FreqShift,MaxCCoef] = MeasureFreqShift( InitMetabRefTimeSerie,RefTimeSerie, Time,FreqRange,FreqPrec );
        RefTimeSerie=RefTimeSerie.*exp(2*pi*1i*Time*FreqShift);

        [FreqMapCorrMetab, CCoefMapMetab,RefSpectrumShortMetab ] = MeasureFreqMapNoSmoothing(RefTimeSerie, permute(formTensorProduct(U_rrrc, V_tc*S),[4,1,2,3]),BrainMask,mrsiReconParams);
   
        %Lipids
 
	[~,Sorig,Vorig] = svd(reshape( SkMask_rrrt.*LipidData_rrrt ,[M*N*Slc NbT ]),0);
	A=8;
	Vlip_tc=Vorig(:,1:A);
	Slip=Sorig(1:A,1:A);
	RefTimeSerie=conj(sum((Vlip_tc*Slip),2));
	if (CorrB0Map_count==1)
       	 InitLipidRefTimeSerie=RefTimeSerie;
	end
	
	Fs = mrsiReconParams.mrProt.samplerate*numel(RefTimeSerie)/mrsiReconParams.mrProt.VSize;
	Time = ([0 :(numel(RefTimeSerie)-1)]'/(Fs));
	[FreqShift,MaxCCoef] = MeasureFreqShift( InitLipidRefTimeSerie,RefTimeSerie, Time,FreqRange,FreqPrec );
	RefTimeSerie=RefTimeSerie.*exp(2*pi*1i*Time*FreqShift);
	
	[FreqMapCorrLip, CCoefMapLip,RefSpectrumShortLip ] = MeasureFreqMapNoSmoothing(RefTimeSerie, permute(LipidData_rrrt,[4,1,2,3]),mrsiReconParams.SkMask,mrsiReconParams);

        FreqMapCorr = FreqMapCorrMetab + FreqMapCorrLip;

        %smoothing of Maps
        
        sigma=mrsiReconParams.GaussianSigma;
        [X,Y,Z] = ndgrid(1:M, 1:N,1:Slc);
        xc=floor(M/2)+1;yc=floor(N/2)+1;zc=floor(Slc/2)+1;
        exponent = -((X-xc).^2 + (Y-yc).^2 + (Z-zc).^2)./(2*sigma^2);
        Kernel = fftshift(fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2),3);
        Kernel = fft(fft(fft(Kernel,[],1),[],2),[],3);
        FreqMapCorr=real(ifft(ifft(ifft(Kernel.*fft(fft(fft(FreqMapCorr,[],1),[],2),[],3),[],1),[],2),[],3));

        MFreqMapCorr=mean(abs(FreqMapCorr(BrainMask>0)));
        MFreqMap=mean(abs(FreqMap(BrainMask>0)));
        RelFreqCorr=MFreqMapCorr/MFreqMap;
        fprintf('Frequency Map correction Nb.%g : RelFreqCorr = %g, MeanFreqCorr= %g, MeanFreqMap= %g, it = %g \n', CorrB0Map_count,RelFreqCorr,MFreqMapCorr,MFreqMap,k+1);

        if(RelFreqCorr>0.01) %If the Freq correction is significant

            FreqMap=single(FreqMap-FreqMapCorr.*mrsiReconParams.ImMask);

            FreqshiftCorr_rrrt=exp(2*pi*1i*Time_111t.*FreqMapCorr.*mrsiReconParams.ImMask);
            MGrad_rrrt = MGrad_rrrt.*FreqshiftCorr_rrrt;
            RGradLipRM_rrrt = RGradLipRM_rrrt.*FreqshiftCorr_rrrt;
            
            UpdateMetabData_rrrt = UpdateMetabData_rrrt.*FreqshiftCorr_rrrt;
    	     MetabData_rrrt = MetabData_rrrt.*FreqshiftCorr_rrrt;
     
            UpdateLipidData_rrrt = UpdateLipidData_rrrt.*FreqshiftCorr_rrrt;
            LipidData_rrrt = LipidData_rrrt.*FreqshiftCorr_rrrt;
            clear FreqshiftCorr_rrrt;
            
            Lipid_stack_rf = fft(reshape(LipidData_rrrt(SkMask_rrrt>0),[],NbT),[],2);

		LipReg=LipReg/10;
		diagLRM=1;
		LipExp=1;
		while diagLRM>TolLipRMDamping
	   		if LipReg*1.2<5E5
			    LipReg=LipReg*1.2;
		    	else
		    	    LipExp=LipExp*1.05;
		    	end
		    LipidRM_ff = (inv( eye(NbT) + LipReg * Lipid_stack_rf'*(Lipid_stack_rf)/norm(Lipid_stack_rf)^2 ))^LipExp;
		    diagLRM=mean(abs(diag(LipidRM_ff)));

		end

			    
            s=[ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_FreShiftCorr_step', num2str(k+1), '.ps'];
            if exist(s,'file');delete(s);end

            figs=figure('visible','off');
            imagesc(Vol2Image(FreqMapCorr));colorbar;
            title('Frequency Map correction');
            print(figs, '-append', '-dpsc2',s);

            imagesc(Vol2Image(FreqMap));colorbar;
            title('Resulting Frequency Map');
            print(figs, '-append', '-dpsc2',s);

            imagesc(Vol2Image(FreqMapCorrMetab));colorbar;
            title('Frequency Map correction from Metabolite signal');
            print(figs, '-append', '-dpsc2',s);


            imagesc(Vol2Image(CCoefMapMetab));colorbar;
            title('Correlation Map from Metabolite signal');
            print(figs, '-append', '-dpsc2',s);
            plot(1:numel(RefSpectrumShortMetab),real(RefSpectrumShortMetab),1:numel(RefSpectrumShortMetab),imag(RefSpectrumShortMetab));
            title('Reference Cropped Spectrum from Metabolite signal');
            print(figs, '-append', '-dpsc2',s);
            

	    figs=figure('visible','off');
	    imagesc(Vol2Image(FreqMapCorrLip));colorbar;
	    title('Frequency Map correction from Lipid signal');
	    print(figs, '-append', '-dpsc2',s);
	    
	    
	    imagesc(Vol2Image(CCoefMapLip));colorbar;
	    title('Correlation Map  from Lipid signal');
	    print(figs, '-append', '-dpsc2',s);
	    plot(1:numel(RefSpectrumShortLip),real(RefSpectrumShortLip),1:numel(RefSpectrumShortLip),imag(RefSpectrumShortLip));
	    title('Reference Cropped Spectrum  from Lipid signal');
	    print(figs, '-append', '-dpsc2',s);


            close
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep other norms for monitoring


    norm_xi=norm(xi_(:));
    norm_u=norm(U_rrrc(:));%norm(ww(:) + divp(:));
    norm_r=norm(NormGrad(:));
    norm_divp=norm(divp(:));
    norm_ww=norm(ww(:));


    StepDiffu = [ StepDiffu norm_u];
    StepDiffr = [ StepDiffr norm_r];
    StepDiffxi = [ StepDiffxi norm_xi];
    StepDiffq = [ StepDiffq norm_q];
    StepDiffp = [ StepDiffp norm_p];
    StepDiffww = [ StepDiffww norm_ww ];
    StepDiffdivp = [ StepDiffdivp norm_divp ];


    if mod(k+1,check_it) == 0

        dataConsistencyCost =  norm(Grad_rrrt(:));
        StepDiff = [ StepDiff, (dataConsistencyCost-costFunVal(end))/(dataConsistencyCost*check_it)];
        costFunVal = [costFunVal, dataConsistencyCost];

        fprintf('\nTGV2-L2-3D: it = %g, costFunDiff = %g,tau_p = %g, RelDiffSq_U = %g, Alpha0 = %g ', k+1,StepDiff(end) ,tau_p, RelDiffSq_U,alpha0);
        fprintf('\nSpectral Iteration: %d , RelDiffSq_V = %d, step_noConv = %d', k+1, RelDiffSq_V,step_noConv);

        CovV=V_tc'*V_tc;
        U_rc =reshape(U_rrrc,[],size(U_rrrc,4));
        CovU=(U_rc'*U_rc);

        fprintf('\nTime component independance: %d , Spatial component independance: %d\n', trace(abs(CovV))/sum(abs(CovV(:))),  trace(abs(CovU))/sum(abs(CovU(:))));    
    
	EnergyLR=(reshape(MetabDataLR_rrrf+ fft(LipidData_rrrt,[],4),[],NbT) );
	EnergyLR = sqrt(sum(abs(EnergyLR(:,FreqOI)).^2,2));
	NoLRNRG = median(EnergyLR(mrsiReconParams.SkMask>0));
        EnergyLR=(reshape(MetabDataLR_rrrf+ fft(LipidData_rrrt,[],4),[],NbT) )*LipidRM_ff*SBRMOp_ff;
	EnergyLR = sqrt(sum(abs(EnergyLR(:,FreqOI)).^2,2));
	LRNRG = median(EnergyLR(mrsiReconParams.SkMask>0)) ;
	Ratio = (NoLRNRG/LRNRG-1);
	
        fprintf(['k=',num2str(k),', NormGrad=',num2str(norm(RGradMetabShort_rrrt(:))),', ']);
    	fprintf(['tau_p=',num2str(tau_p),', Ratio=',num2str(Ratio),', LipReg=',num2str(LipReg),', LipExp=',num2str(LipExp),', RelDiffMetab=',num2str(RelDiffMetab),', RelDiffLipid=',num2str(RelDiffLipid),'\n']); 

    end 
    if  mod(k+1,Plot_it)==0

        fprintf('Make plots ...');
	[Uorig,Sorig,Vorig] = svd(reshape( SkMask_rrrt.*LipidData_rrrt ,[M*N*Slc NbT ]),0);
	Vlip_tc=Vorig(:,1:Nlip);
	Slip=Sorig(1:Nlip,1:Nlip);
	Ulip_rrrc=Uorig(:,1:Nlip);
	Ulip_rrrc = SignalMask_rrr1.*reshape(Ulip_rrrc,[M N Slc Nlip]);
	
        VisualizeTGV( Ulip_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_step', num2str(k+1),'_Lipid']);
        VisualizeSpectral( Vlip_tc,Slip, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_step', num2str(k+1),'_Lipid'])
        VisualizeTGV( U_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_step', num2str(k+1),'_Metabs']);
        VisualizeSpectral( V_tc,S, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_step', num2str(k+1),'_Metabs'])
      
    end

end

%Gather memory from the GPU
if thisgpu>0
	VarNames=who;
	for VN=1:numel(VarNames)
		 eval([VarNames{VN}, ' = gatherFromGPU(',VarNames{VN},');']);
	end
end

%Reorder Component following signular value:
[SDiag, DescOrder]=sort(diag(S),'descend');
U_rrrc=U_rrrc(:,:,:,DescOrder).*repmat(mrsiReconParams.ImMask,[1 1 1 NbComp]);
V_tc=V_tc(:,DescOrder);
S=diag(SDiag);

fprintf([ '\nTGV Recon done in ', num2str(k),' steps. Relative Step Diff U =',  num2str(RelDiffSq_U), ' steps. Relative Step Diff V =',  num2str(RelDiffSq_V),'\n']);
