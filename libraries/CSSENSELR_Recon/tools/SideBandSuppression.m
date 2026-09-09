function [Data_fkkk,FreqShift] = SideBandSuppression(Data_fkkk,PtRange,Rank,mrsiReconParams)
NbPt=(PtRange(2)-PtRange(1)+1);
DownField_tk = ifft(Data_fkkk((end-PtRange(2)+1):(end-PtRange(1)+1),:,:,:),[],1);
DownField_tk = reshape(DownField_tk,NbPt,[]);
[~,~,V] = svd(permute(DownField_tk,[2,1]),0);
V_tc=conj(V(:,1:Rank));
SideBands_tk=(V_tc*V_tc')*DownField_tk;

clear DownField_tk

%Margin=round(NbPt/4);
UpField_tk =  ifft(Data_fkkk(PtRange(1):PtRange(2),:,:,:),[],1);
UpField_tk = reshape(UpField_tk,NbPt,[]);
Fs=mrsiReconParams.mrProt.samplerate*NbPt/mrsiReconParams.mrProt.VSize;
Time=[0 :(NbPt-1)]'/Fs;
FreqRange=40;,FreqPrec=1;
CleanDataUpField_tk = 0*UpField_tk;
FreqShift=zeros(1,size(SideBands_tk,2));
IOP=eye(NbPt);
VShftd_tc=0*V_tc;
 fprintf([ 'Processing point']); 
for a=1:size(SideBands_tk,2)
	if(mod(a,round(size(SideBands_tk,2)/100))==0)
 		fprintf([ ' ', num2str(round(100*a/size(SideBands_tk,2))), '%%,']);drawnow('update');
    end	   
   [FreqShift(a),MaxCCoef] = MeasureFreqShift( conj(SideBands_tk(:,a)),UpField_tk(:,a), Time,FreqRange,FreqPrec );  
   VShftd_tc=conj(V_tc).*exp(-2*pi*1i*Time*FreqShift(a));
   CleanDataUpField_tk(:,a)=(IOP -VShftd_tc*VShftd_tc')*UpField_tk(:,a);
   
end
fprintf('\n');

% VShftd_tc=conj(V_tc);
% CleanDataUpField_tk=(eye(NbPt)-VShftd_tc*VShftd_tc')*UpField_tk;

clear VShftd_tc UpField_tk SideBands_tk


FreqShift = reshape(FreqShift,[size(Data_fkkk,2),size(Data_fkkk,3),size(Data_fkkk,4)]);
Data_fkkk((PtRange(1):PtRange(2)),:,:,:) = reshape(fft(CleanDataUpField_tk,[],1),[NbPt,size(Data_fkkk,2),size(Data_fkkk,3),size(Data_fkkk,4)]);
end

