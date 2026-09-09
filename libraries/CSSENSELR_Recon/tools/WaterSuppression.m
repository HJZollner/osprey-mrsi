function [filtMRSI ]=WaterSuppression(MRSIData_tkkk,mrProt,FiltParam,kmask)


filtMRSI=zeros(size(MRSIData_tkkk));
step=0;
tot_steps= sum(kmask(:));
for c=1:size(MRSIData_tkkk,4);
    for b=1:size(MRSIData_tkkk,3);
        for a=1:size(MRSIData_tkkk,2);
            if(kmask(a,b,c))
                step=step+1;
                %             if(mod(step,round(tot_steps/5))==0)
                %                 fprintf(sprintf('filter step %g out of %g. \n',step,tot_steps));
                %             end
                filtMRSI(:,a,b,c) =  Fast_HSVD_Filter(MRSIData_tkkk(:,a,b,c),1.0/mrProt.DwellTime,FiltParam.Comp,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
            end
        end
    end
end

end