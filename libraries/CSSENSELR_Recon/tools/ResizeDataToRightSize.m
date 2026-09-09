function Resized_raw_tckkk=ResizeDataToRightSize(raw_tckkk,NPhase,NRead,NSlice);

[NbT NbC M N S]=size(raw_tckkk);
data_kkk=squeeze(sum(sum(abs(raw_tckkk),1),2));

[ maxValue, Imax]=max(data_kkk(:));
[CM_r CM_p CM_s ] =ind2sub(size(data_kkk),Imax);

P_CM_r=log2(CM_r-1);
P_CM_p=log2(CM_p-1);
P_CM_s=log2(CM_s-1);
norm_conv_cm=abs(P_CM_r-round(P_CM_r))/(P_CM_r) +abs(P_CM_p-round(P_CM_p))/(P_CM_p)+abs(P_CM_s-round(P_CM_s))/(P_CM_s);
if(norm_conv_cm > 1e-4)
	warning('Original k-space center of mass (2^n+1 x 2^m +1 x 2^p +1 expected) is unusual in ResizeDataToRightSize.m!');
        fprintf(['Original center of mass is:[',num2str(CM_r),' ',num2str(CM_p),' ',num2str(CM_s),'].\n']);
        CM_r=2^round(P_CM_r)+1; CM_p=2^round(P_CM_p)+1;CM_s=2^round(P_CM_s)+1;
        fprintf(['Center of mass was corrected to:[',num2str(CM_r),' ',num2str(CM_p),' ',num2str(CM_s),'].\n']);  
end

norm_k1=sum(sum(abs(data_kkk),2),3);
min_k1=min(find([norm_k1>0]));
max_k1=max(find([norm_k1>0]));

norm_k2=sum(sum(abs(data_kkk),1),3);
min_k2=min(find([norm_k2>0]));
max_k2=max(find([norm_k2>0]));

norm_k3=sum(sum(abs(data_kkk),1),2);
min_k3=min(find([norm_k3>0]));
max_k3=max(find([norm_k3>0]));

raw_tckkk=raw_tckkk(:,:,min_k1:max_k1,min_k2:max_k2,min_k3:max_k3);
[NbT NbC M N S]=size(raw_tckkk);

CM_r=CM_r-min_k1+1;
CM_p=CM_p-min_k2+1;
CM_s=CM_s-min_k3+1;


GoalCM_r=round((NRead +1)/2);
GoalCM_p=round((NPhase+1)/2);
GoalCM_s=round((NSlice+1)/2);

shift_r=(GoalCM_r-CM_r);
shift_p=(GoalCM_p-CM_p);
shift_s=(GoalCM_s-CM_s);

Resized_raw_tckkk=zeros([NbT NbC NRead NPhase NSlice]);

Resized_raw_tckkk(:,:,(1+shift_r):(M+shift_r),(1+shift_p):(N+shift_p),(1+shift_s):(S+shift_s))=raw_tckkk;

Resized_raw_tckkk=permute(Resized_raw_tckkk,[1 2 4 5 3]);

end
