function [U_rrrc,V_tc, S,costFunVal,FreqMap ] = tgv2_l2_3D_SENSE_MRSILowRank(Data_ckkkt,U_rrrc,V_tc,S,OriginalData_frrr, alpha0, alpha1,maxits,minits,mrsiReconParams,Threshold)
% Primal dual TGV2 algorithm, as described in the TGV paper

Data_ckkkt=(single(Data_ckkkt));
U_rrrc=(single(U_rrrc));
V_tc=(single(V_tc));
S=(single(S));

[dxm,dym,dzm,dxp,dyp,dzp] = defineDiffOperators();

check_it =mrsiReconParams.LRTGVModelParams.check_it;
Plot_it=mrsiReconParams.LRTGVModelParams.Plot_it;
CorrB0Map_it=mrsiReconParams.LRTGVModelParams.CorrB0Map_it;
CorrB0Map_Maxcount=mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount;
Orthogonalize_it=mrsiReconParams.LRTGVModelParams.Orthogonalize_it;
reduction =mrsiReconParams.LRTGVModelParams.reduction;
min_SpectStep=mrsiReconParams.LRTGVModelParams.min_SpectStep;
max_SpectStep=mrsiReconParams.LRTGVModelParams.max_SpectStep;
min_taup=mrsiReconParams.LRTGVModelParams.min_taup;
max_taup=mrsiReconParams.LRTGVModelParams.max_taup;
DualPrimalTauFact=mrsiReconParams.LRTGVModelParams.DualPrimalTauFact;
CorrB0Map_count=1;

Init_U_rrrc= U_rrrc;
[ M N Slc NbComp] = size(U_rrrc); % numSamplesOnSpoke, numSamplesOnSpoke, nCh
NbCoil = size(Data_ckkkt,1);
NbT= size(Data_ckkkt,5);
SizeVol = size(U_rrrc);
DimVol = ndims(U_rrrc)-1;
UIndcs  = repmat({':'}, [1, numel(SizeVol)]);
numSpatialPts=size(U_rrrc, 1)*size(U_rrrc,2)*size(U_rrrc, 3);


%for Forward Transform
FreqMap=mrsiReconParams.WaterFreqMap ;
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
DelayT=mrsiReconParams.AcqDelay;
Time_rrrt=single(permute(repmat(([0 :(NbT-1)]'/Fs+DelayT),[1 M N Slc]),[2,3,4,1]));
%Freqshift_rrrt=exp(2*pi*1i*Time_rrrt.*repmat(FreqMap,[1 1 1 NbT]));
Freqshift_1rrrt=(single(reshape(exp(2*pi*1i*Time_rrrt.*repmat(FreqMap,[1 1 1 NbT])),[1 M N Slc NbT])));
SENSE_crrr1=(single(reshape(mrsiReconParams.SENSE,[size(mrsiReconParams.SENSE) 1])));
kmask_1kkk1=(single(reshape(mrsiReconParams.kmask,[1 size(mrsiReconParams.kmask) 1])));
HannF_1kkk1=(single(reshape(mrsiReconParams.HKernel,[1 size(mrsiReconParams.HKernel) 1])));
BMask_rrr1=(single(reshape(mrsiReconParams.BrainMask,[size(mrsiReconParams.BrainMask) 1])));
SENSE_crrr1=(single(SENSE_crrr1.* reshape(BMask_rrr1,[1 size(BMask_rrr1)])));

G=0*Freqshift_1rrrt;

StepNormGrad = [];
StepDiff=Threshold;
StepDiffu = [];
StepDiffr = [];
StepDiffxi = [];
StepDiffq = [ ];
StepDiffp = [];
StepDiffww = [ ];
StepDiffdivp = [];
StepTauP = [];
StepStepSize = [];

tau_p = max_taup;
tau_d= tau_p*2;
stepSize=max_SpectStep;


PrevNormGradU=1E12;
NormGradWentDown=0;
NbDivergSteps=0;

RelDiffSq_V=1;
RelDiffSq_U=1;

t_old    = 2;
PrevNormGrad=0;
step_noConv=0;

p = single(zeros([SizeVol,3]));
q = single(zeros([SizeVol,6]));
xi = single(zeros([SizeVol,3]));

U_rrrc=U_rrrc.*BMask_rrr1;
u = U_rrrc;
xi_ = xi; % v in article
V_old    = V_tc;
V_proj    = V_tc;
uold = u; xiold = xi;


G  = reshape(formTensorProduct(U_rrrc, V_tc*S),[1 M N Slc NbT]); %1rrrt
DataError = fft(fft(fft(Freqshift_1rrrt.*SENSE_crrr1.*G,[],2),[],3),[],4).*HannF_1kkk1.*kmask_1kkk1; %ckkkt
dataConsistencyCost = norm(Data_ckkkt(:) - DataError(:) );
DataFidelCostVal = dataConsistencyCost;

TGVCost0Val = 0;
TGVCost1Val = 0;


k=-1;

alpha00 = alpha0/reduction;
alpha10 = alpha1/reduction;
alpha002 = alpha0*reduction;
alpha102 = alpha1*reduction;
alpha01 = alpha0;
alpha11 = alpha1;

%store variable to GPU
Vars=whos;
thisgpu =SelectFreeGPU(1.5*sum([Vars.bytes]));

if thisgpu>0
	fprintf("Running Recon on GPU nb: "+num2str(thisgpu)+" ...\n");
	VarNames=who;
	for VN=1:numel(VarNames)
		 eval([VarNames{VN}, ' = transferToGPU(',VarNames{VN},');']);
    end
end

%tic
while k<maxits
    k=k+1;
    %toc
    %tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL CONVERGENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % update alpha's

    if k<=(minits/2)
     	 alpha0 = exp(k/(minits/2)*log(alpha01) + ((minits/2)-k)/(minits/2)*log(alpha00));
      	 alpha1 = exp(k/(minits/2)*log(alpha11) + ((minits/2)-k)/(minits/2)*log(alpha10));
    elseif k<=(minits)
        alpha0 = exp((k-minits/2)/(minits/2)*log(alpha01) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha002));
     	 alpha1 = exp((k-minits/2)/(minits/2)*log(alpha11) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha102));
    else
   	 alpha0 = alpha01;
     alpha1 = alpha11;
    end

    if k > minits
      	k_ind=(k-minits)/(maxits-minits); %from 0 to 1
    	max_SpectStep=mrsiReconParams.LRTGVModelParams.max_SpectStep;
    	max_taup=mrsiReconParams.LRTGVModelParams.max_taup;

    	max_taup = min_taup + (max_taup-min_taup)*exp(-k_ind*10);
    	max_SpectStep = min_SpectStep +  (max_SpectStep-min_SpectStep)*exp(-k_ind*10);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uold = u;
    xiold = xi;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DUAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % data fidelity
    G  = reshape(formTensorProduct(U_rrrc, V_tc*S),[1 M N Slc NbT]); %1rrrt
    DataError = (fft(fft(fft(Freqshift_1rrrt.*SENSE_crrr1.*G,[],2),[],3),[],4).*HannF_1kkk1 - Data_ckkkt).*kmask_1kkk1;%c-k-k-k-t
    GeneralGrad = BMask_rrr1.*squeeze(sum(conj(SENSE_crrr1).*conj(Freqshift_1rrrt).*ifft(ifft(ifft( DataError,[],2),[],3),[],4),1));% r-r-r-t


    % gradient
    ux = dxp(U_rrrc);
    uy = dyp(U_rrrc);
    uz = dzp(U_rrrc);

    p(UIndcs{:},1) = p(UIndcs{:},1) - tau_d*(ux + xi_(UIndcs{:},1));
    p(UIndcs{:},2) = p(UIndcs{:},2) - tau_d*(uy + xi_(UIndcs{:},2));
    p(UIndcs{:},3) = p(UIndcs{:},3) - tau_d*(uz + xi_(UIndcs{:},3));

    % projection

    for comp=1:NbComp
        absp = sqrt(abs(p(UIndcs{1:(end-1)},comp,1)).^2 + abs(p(UIndcs{1:(end-1)},comp,2)).^2  + abs(p(UIndcs{1:(end-1)},comp,3)).^2);
        denom = max(1,absp/(alpha1*max(diag(S))/S(comp,comp)));
        p(UIndcs{1:(end-1)},comp,1) = p(UIndcs{1:(end-1)},comp,1)./denom;
        p(UIndcs{1:(end-1)},comp,2) = p(UIndcs{1:(end-1)},comp,2)./denom;
        p(UIndcs{1:(end-1)},comp,3) = p(UIndcs{1:(end-1)},comp,3)./denom;
    end

    % symmetrized gradient
    g_xi1 = dxm(xi_(UIndcs{:},1));
    g_xi2 = dym(xi_(UIndcs{:},2));
    g_xi3 = dzm(xi_(UIndcs{:},3));
    g_xi4 = (dxm(xi_(UIndcs{:},2)) + dym(xi_(UIndcs{:},1)))/2;
    g_xi5 = (dxm(xi_(UIndcs{:},3)) + dzm(xi_(UIndcs{:},1)))/2;
    g_xi6 = (dym(xi_(UIndcs{:},3)) + dzm(xi_(UIndcs{:},2)))/2;

    q(UIndcs{:},1) = q(UIndcs{:},1) - tau_d*g_xi1; % line
    q(UIndcs{:},2) = q(UIndcs{:},2) - tau_d*g_xi2;
    q(UIndcs{:},3) = q(UIndcs{:},3) - tau_d*g_xi3;
    q(UIndcs{:},4) = q(UIndcs{:},4) - tau_d*g_xi4; % line
    q(UIndcs{:},5) = q(UIndcs{:},5) - tau_d*g_xi5;
    q(UIndcs{:},6) = q(UIndcs{:},6) - tau_d*g_xi6;

    % projection

    for comp=1:NbComp
        absq = sqrt(abs(q(UIndcs{1:(end-1)},comp,1)).^2 + abs(q(UIndcs{1:(end-1)},comp,2)).^2 + abs(q(UIndcs{1:(end-1)},comp,3)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,4)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,5)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,6)).^2);
        denom = max(1,absq/(alpha0*max(diag(S))/S(comp,comp)));
        q(UIndcs{1:(end-1)},comp,1) = q(UIndcs{1:(end-1)},comp,1)./denom;
        q(UIndcs{1:(end-1)},comp,2) = q(UIndcs{1:(end-1)},comp,2)./denom;
        q(UIndcs{1:(end-1)},comp,3) = q(UIndcs{1:(end-1)},comp,3)./denom;
        q(UIndcs{1:(end-1)},comp,4) = q(UIndcs{1:(end-1)},comp,4)./denom;
        q(UIndcs{1:(end-1)},comp,5) = q(UIndcs{1:(end-1)},comp,5)./denom;
        q(UIndcs{1:(end-1)},comp,6) = q(UIndcs{1:(end-1)},comp,6)./denom;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIMAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dual operator

    ww = formTensorProduct(GeneralGrad ,(V_tc*inv(S))' ) ;%r-r-r-c

    % divergence
    divp = dxm(p(UIndcs{:},1)) + dym(p(UIndcs{:},2)) + dzm(p(UIndcs{:},3));

    u = BMask_rrr1.*(u - tau_p*(ww + divp)); %lines 8

    % divergence
    divq1 = dxp(q(UIndcs{:},1)) + dyp(q(UIndcs{:},4)) + dzp(q(UIndcs{:},5));
    divq2 = dxp(q(UIndcs{:},4)) + dyp(q(UIndcs{:},2)) + dzp(q(UIndcs{:},6));
    divq3 = dxp(q(UIndcs{:},5)) + dyp(q(UIndcs{:},6)) + dzp(q(UIndcs{:},3));

    xi(UIndcs{:},1) = xi(UIndcs{:},1) - tau_p*(divq1 - p(UIndcs{:},1));%line 11
    xi(UIndcs{:},2) = xi(UIndcs{:},2) - tau_p*(divq2 - p(UIndcs{:},2));%line 11
    xi(UIndcs{:},3) = xi(UIndcs{:},3) - tau_p*(divq3 - p(UIndcs{:},3));%line 11

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIARY UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_rrrc =2*u - uold;
    xi_ = 2*xi - xiold;

    if  mod(k+1,Orthogonalize_it)==0
        U_rrrc = reshape(OrthogonalizeComponents(reshape(U_rrrc,[size(U_rrrc,1)*size(U_rrrc,2)*size(U_rrrc,3) , size(U_rrrc,4)])),size(U_rrrc));
        u = reshape(OrthogonalizeComponents(reshape(u,[size(u,1)*size(u,2)*size(u,3) , size(u,4)])),size(u));
    end

    for c=1:NbComp
        NormUrrc=sqrt(sum(sum(sum(abs(U_rrrc(:,:,:,c)).^2,1),2),3));
        S(c,c)=S(c,c)*NormUrrc;
        U_rrrc(:,:,:,c) = U_rrrc(:,:,:,c)/NormUrrc;
        u(:,:,:,c) = u(:,:,:,c)/NormUrrc;
    end


    RelDiffSq_U = norm(U_rrrc(:)-u(:))^2/norm(u(:))^2 + norm(xi_(:)-xi(:))^2/norm(xi(:))^2;
    norm_p=norm(p(:),1);
    norm_q=norm(q(:),1);
    RelDiffSq_old=RelDiffSq_U;

    NormGradU=norm(ww(:) + divp(:));


    if k>10 & NormGradU<PrevNormGradU;
        NormGradWentDown=1;
        tau_p = tau_p*1.04;%1.1
        if tau_p>max_taup; tau_p=max_taup;end;
    elseif NormGradWentDown==1;
        tau_p=tau_p*0.95;%0.85
        if tau_p<min_taup;tau_p=min_taup;end;
    end
    tau_d=tau_p*DualPrimalTauFact;
    PrevNormGradU=NormGradU;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECTRAL CONVERGENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    gradStep = ((reshape(U_rrrc, [], NbComp)*inv(S))' * reshape(GeneralGrad , [numSpatialPts, NbT]))' ;

    [ gradStep] = OrthogonalizeComponents(gradStep);



    V_proj   = V_tc -stepSize * gradStep;

    for c=1:size(V_tc,2)
        %if (norm(V_proj(:,c)) > 1)
        V_proj(:,c) = V_proj(:,c) ./ norm(V_proj(:,c));
        %end
    end

    if k>0
        RelDiffSq_V=norm(V_proj(:)-V_old(:))^2/norm(V_proj(:))^2;
    end
    NormGrad=norm(gradStep(:));
    if (NormGrad<PrevNormGrad || k==0)
        step_noConv=0;
        stepSize = 1.03 * stepSize;
        if stepSize>max_SpectStep;stepSize=max_SpectStep;end
    else
        step_noConv=step_noConv+1;
        stepSize = 0.95 * stepSize;
        if stepSize<min_SpectStep;stepSize=min_SpectStep;end
    end

    PrevNormGrad=NormGrad;
    StepNormGrad = [StepNormGrad NormGrad];

    t_new    = (1 + sqrt(1 + 4*t_old.^2)) / 2;
    fact=((t_old - 1) / t_new);

    V_tc = (1-fact)*V_old +(fact)* V_proj;
    for c=1:size(V_tc,2)
        V_tc(:,c) = V_tc(:,c) ./ norm(V_tc(:,c));
    end

    V_old    = V_proj;
    t_old    = t_new;


    if  mod(k+1,Orthogonalize_it)==0
        V_tc = OrthogonalizeComponents(V_tc);
        V_old= OrthogonalizeComponents(V_old);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FREQUENCY MAP DYNAMIC CORRECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if mod(k+1,CorrB0Map_it) == 5 & CorrB0Map_count<=CorrB0Map_Maxcount %& k>100
    	CorrB0Map_count=CorrB0Map_count+1;

        A=8;
        RefTimeSerie=conj(sum((V_tc(:,1:A)*S(1:A,1:A)),2));
        [FreqMapCorr, CCoefMap,RefSpectrumShort ] = MeasureFreqMap(RefTimeSerie, permute(formTensorProduct(U_rrrc, V_tc*S),[4,1,2,3]),mrsiReconParams);
        MFreqMapCorr=mean(abs(FreqMapCorr(mrsiReconParams.BrainMask>0)));
        MFreqMap=mean(abs(FreqMap(mrsiReconParams.BrainMask>0)));
        RelFreqCorr=MFreqMapCorr/MFreqMap;
        fprintf('Frequency Map correction Nb.%g : RelFreqCorr = %g, MeanFreqCorr= %g, MeanFreqMap= %g, it = %g \n', CorrB0Map_count-1,RelFreqCorr,MFreqMapCorr,MFreqMap,k+1);

        if(RelFreqCorr>0.01)

            FreqMap=FreqMap-FreqMapCorr.*mrsiReconParams.BrainMask;

    	    Freqshift_1rrrt=reshape(exp(2*pi*1i*Time_rrrt.*repmat(FreqMap,[1 1 1 NbT])),[1 M N Slc NbT]);
            s=[ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_FreShiftCorr_step', num2str(k+1), '.ps'];
            if exist(s);delete(s);end
            figs=figure('visible','off');
            imagesc(Vol2Image(FreqMapCorr));colorbar;
            title('Frequency Map correction');
            print(figs, '-append', '-dpsc2',s);
            imagesc(Vol2Image(FreqMap),[-50 50]);colorbar;
            title('Resulting Frequency Map');
            print(figs, '-append', '-dpsc2',s);

            imagesc(Vol2Image(CCoefMap));colorbar;
            title('Correlation Map');
            print(figs, '-append', '-dpsc2',s);
            plot(1:numel(RefSpectrumShort),real(RefSpectrumShort),1:numel(RefSpectrumShort),imag(RefSpectrumShort));
            title('Reference Cropped Spectrum');
            print(figs, '-append', '-dpsc2',s);

            close
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep other metrics for monitoring


    norm_xi=norm(xi_(:));
    norm_u=norm(U_rrrc(:));
    norm_r=norm(NormGradU(:));
    norm_divp=norm(divp(:));
    norm_ww=norm(ww(:));



    StepDiffu = [ StepDiffu norm_u];
    StepDiffr = [ StepDiffr norm_r];
    StepDiffxi = [ StepDiffxi norm_xi];
    StepDiffq = [ StepDiffq norm_q];
    StepDiffp = [ StepDiffp norm_p];
    StepDiffww = [ StepDiffww norm_ww ];
    StepDiffdivp = [ StepDiffdivp norm_divp ];
    StepTauP = [StepTauP tau_p];
    StepStepSize = [StepStepSize stepSize];

    if mod(k+1,check_it) == 0

    	G  = reshape(formTensorProduct(U_rrrc, V_tc*S),[1 M N Slc NbT]); %1rrrt
        DataError = fft(fft(fft(Freqshift_1rrrt.*SENSE_crrr1.*G,[],2),[],3),[],4).*HannF_1kkk1.*kmask_1kkk1; %ckkkt
    	dataConsistencyCost = norm(Data_ckkkt(:) - DataError(:) );

        StepDiff = [ StepDiff, (dataConsistencyCost-DataFidelCostVal(end))/(dataConsistencyCost*check_it)];
        DataFidelCostVal = [DataFidelCostVal, dataConsistencyCost];

        SNorm_Op_111c1 = reshape(diag(S),[1 1 1 NbComp 1])/max(diag(S));
        tgv1Cost           = sum(abs(vectorizeArray((cat(5, ux, uy, uz) - xi).*SNorm_Op_111c1)));
        tgv0Cost	 = sum(vectorizeArray((abs(g_xi1) + abs(g_xi2) + abs(g_xi3) + 2*abs(g_xi4) + 2*abs(g_xi5) + 2*abs(g_xi6)).*SNorm_Op_111c1));

        TGVCost0Val      = [ TGVCost0Val, tgv0Cost];  % TGV penality term
        TGVCost1Val      = [ TGVCost1Val, alpha1/alpha0*tgv1Cost];  % TGV penality term

        fprintf('\nTGV2-L2-3D: it = %g, costFunDiff = %g,tau_p = %g, RelDiffSq_U = %g, Alpha0 = %g ', k+1,StepDiff(end) ,tau_p, RelDiffSq_U,alpha0);
        fprintf('\nSpectral Iteration: %d , step_size = %d, RelDiffSq_V = %d, step_noConv = %d', k+1, stepSize,RelDiffSq_V,step_noConv);
        CovV=V_tc'*V_tc;
        U_rc =reshape(U_rrrc,[],size(U_rrrc,4));
        CovU=(U_rc'*U_rc);
        fprintf('\nTime component independance: %d , Spatial component independance: %d\n', trace(abs(CovV))/sum(abs(CovV(:))),  trace(abs(CovU))/sum(abs(CovU(:))));


    end
    if  mod(k+1,Plot_it)==0

        subplot(3,3,1);plot(DataFidelCostVal);title('Data Fidelity cost');
        subplot(3,3,2);plot(TGVCost0Val);title('TGV 0 term cost');
        subplot(3,3,3);plot(TGVCost1Val);title('TGV 1 term cost');
        subplot(3,3,4);plot(StepDiffww);title('Data Fid. Spat. Gradient norm');
        subplot(3,3,5);plot(StepNormGrad);title('Data Fid. Spect. Gradient norm');
        subplot(3,3,6);plot(StepTauP);title('Spatial step size');
        subplot(3,3,7);plot(StepStepSize);title('Spectral step size');


        s=[ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_CostValues_Diagnostic.ps'];
        if exist(s);delete(s);end
        print(s, '-dpsc2');

        VisualizeTGV( U_rrrc,[mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_step', num2str(k+1)]);

        VisualizeSpectral( V_tc,S, [ mrsiReconParams.Log_Dir,filesep,'LowRankTGV_Recon',filesep,mrsiReconParams.NameData,'_step', num2str(k+1)])

    end

end

%Compute final cost functions
clear costFunVal
G  = reshape(formTensorProduct(U_rrrc, V_tc*S),[1 M N Slc NbT]); %1rrrt
DataError = fft(fft(fft(Freqshift_1rrrt.*SENSE_crrr1.*G,[],2),[],3),[],4).*HannF_1kkk1.*kmask_1kkk1; %ckkkt
costFunVal(1)  = norm((Data_ckkkt(:) - DataError(:))); % Data Fidelity


SNorm_Op_111c1 = reshape(diag(S),[1 1 1 NbComp 1]);;
tgv1Cost     = sum(abs(vectorizeArray((cat(5, ux, uy, uz) + xi).*SNorm_Op_111c1)));
tgv0Cost	 = sum(vectorizeArray((abs(g_xi1) + abs(g_xi2) + abs(g_xi3) + abs(g_xi4) + abs(g_xi5) + abs(g_xi6)).*SNorm_Op_111c1));

costFunVal(2)       = alpha0/alpha1 * tgv1Cost +  tgv0Cost;  % First and Second order driv. U cos

%Gather memory from the GPU
if thisgpu>0
	VarNames=who;
	for VN=1:numel(VarNames)
		 eval([VarNames{VN}, ' = gatherFromGPU(',VarNames{VN},');']);
	end
end

%Reorder Component following signular value:
[SDiag, DescOrder]=sort(diag(S),'descend');
U_rrrc=U_rrrc(:,:,:,DescOrder).*repmat(mrsiReconParams.BrainMask,[1 1 1 NbComp]);
V_tc=V_tc(:,DescOrder);
S=diag(SDiag);

fprintf([ '\nTGV Recon done in ', num2str(k),' steps. Data Fidelity cost =',  num2str(costFunVal(1) ), '. TGV penality cost =',  num2str(costFunVal(2) ),'\n']);

