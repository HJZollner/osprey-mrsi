function [MRSIdata_t CCoef Q]=WSVD(MRSIdata_tc, Range_noise);

MRSIdata_t=zeros(size(MRSIdata_tc,1),1);
            
            Noise_tc=MRSIdata_tc(Range_noise(1):Range_noise(2),:);
            Noise_tc=detrend(Noise_tc,'constant');% remove the mean (normally zero but in case of bad baseline)
        
            for c=1:size(Noise_tc,2)for d=1:size(Noise_tc,2)
             Noise_corr(c,d)=mean(squeeze(Noise_tc(:,c).*conj(Noise_tc(:,d))));
               end;end
 
            [X,D] = eig(Noise_corr);% following the notation of Rodgers et al. 2016
            W=D^(-0.5)*X';
         
           % filtMRSI = transpose(mrsiExpFilter(transpose(MRSIdata_tc),0.5,250));
           filtMRSI = MRSIdata_tc;
            for t=1:size(MRSIdata_tc,1)
                S(t,:)=W*transpose(squeeze(filtMRSI(t,:)));
            end
            [V,Sig,U]=svd(S);
         
            Q=Sig(1,1)*V(:,1);
           
            CCoef=inv(W)*U(:,1);
            
            for c=1:size(CCoef,1)
                MRSIdata_t=MRSIdata_t+squeeze(CCoef(c)*MRSIdata_tc(:,c));
            end
        

