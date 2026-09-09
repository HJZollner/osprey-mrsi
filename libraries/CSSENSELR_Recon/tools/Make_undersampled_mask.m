function [US_MASK]=Make_undersampled_mask(expo,Fil_Fact,kmask,CSRelHardRadius)

CompSensingAcc=1./Fil_Fact;
kmask=fftshift(fftshift(fftshift(kmask,1),2),3);
Size_grid=size(kmask);
expo=0.5;
dm1 = Size_grid(1);
dm2 = Size_grid(2);
dm3 = Size_grid(3);

lowd1 = -floor(dm1/2);
uppd1 = round(dm1/2)-1;%(dm1%2) ? dm1/2 : dm1/2-1;
lowd2 = -floor(dm2/2);
uppd2 = round(dm2/2)-1;%(dm2%2) ? dm2/2 : dm2/2-1;
lowd3 = -floor(dm3/2);
uppd3 = round(dm3/2)-1;%(dm3%2) ? dm3/2 : dm3/2-1;

addi_voxel_rad=0.125;
TotProb=0;
d=0;

ia3=1;ia2=1;ia1=1;
CSMask_Filling=0;
Elliptic_Filling=0;
HardRadius_Filling=0;

for a3=lowd3:uppd3;
    ia2=1;
    for a2=lowd2:uppd2;
        ia1=1;
        for a1=lowd1:uppd1
            if(uppd1 == 0) d=0;
            else d= a1;
            end
            
            if (uppd1 >0)
                dist_k = d*d/((uppd1+1)*(uppd1+1));
            end
            
            d = d/((uppd1+addi_voxel_rad));
            dist = d*d;
            
            if(uppd2 == 0) d=0;
            else d=a2;
            end
            if (uppd2 >0)
                dist_k = dist_k + d*d/((uppd2+1)*(uppd2+1));
            end
            d = d/((uppd2+addi_voxel_rad));
            dist = dist + d*d;
            
            if(uppd3 == 0) d=0 ;
            else d= a3;
            end
            
            if (uppd3 >0)
                dist_k = dist_k + d*d/((uppd3+1)*(uppd3+1));
            end
            d = d/((uppd3+addi_voxel_rad));
            dist = dist + d*d;
            
            dist_k =  dist_k -CSRelHardRadius*CSRelHardRadius;
            
            dist = sqrt( dist );
            Mask_DistRel(ia1,ia2,ia3)=dist;
            %compute probbability map with hard center
            ProbMap(ia1,ia2,ia3)=0;
            %if(dist_k<(double)m_iCSHardRadius){
            if(dist_k<0)
                ProbMap(ia1,ia2,ia3)=1000;
                HardRadius_Filling=HardRadius_Filling+1;        
            else
                ProbMap(ia1,ia2,ia3)=1.0-sqrt(dist_k);
                %ProbMap(ia1,ia2,ia3)=1.0/(sqrt(dist_k).^expo);
            end
            if (ProbMap(ia1,ia2,ia3)<0)ProbMap(ia1,ia2,ia3)=0;end
            
            %compute elliptic k-space mask and filling
            Elliptic_Mask(ia1,ia2,ia3)=0;
            %if (dist<=1 )
            if (kmask(ia1,ia2,ia3)==1)
                Elliptic_Mask(ia1,ia2,ia3)=1;
                Elliptic_Filling=Elliptic_Filling+1;
                if (ProbMap(ia1,ia2,ia3)<1)
                    TotProb=TotProb+ProbMap(ia1,ia2,ia3);
                else
                    TotProb=TotProb+1;
                end
            end
            
            ia1=ia1+1;
        end
        ia2=ia2+1;
    end
    ia3=ia3+1;
end

RelFill=TotProb/Elliptic_Filling;
GoalFill=Elliptic_Filling*1.0/CompSensingAcc;

if((Elliptic_Filling*1.0/HardRadius_Filling)<CompSensingAcc)
    error("Error in CS: Acc. factor too high for this hard radius!" );
end

%{
		 std::cout<<"RelFill:" <<RelFill<<std::endl;
		std::cout<<"GoalFill:" <<GoalFill<<std::endl;
		std::cout<<"CompSensingAcc:" <<CompSensingAcc<<std::endl;
		std::cout<<"CSMask_Filling before:" <<CSMask_Filling<<std::endl;
		std::cout<<"Elliptic_Filling:" <<Elliptic_Filling<<std::endl;
		CSMask_Filling = 0;
		int stepCS=0;
		double CorrFact=1.0;
%}
CSMask_Filling = Elliptic_Filling;
CS_Mask=Elliptic_Mask;
stepCS=0;
CorrFact=1.0;
while(  abs(GoalFill-CSMask_Filling)*1.0/GoalFill>0.01 && stepCS<100) % precise up to 1%
    %while(abs(GoalFill-CSMask_Filling)*1.0/GoalFill>0.01 ){ /// precise up to 1%
    
    stepCS=stepCS+1;
    CSMask_Filling = Elliptic_Filling;



    ia3=1;
    for a3=lowd3:uppd3
        ia2=1;
        for a2=lowd2:uppd2
            ia1=1;
            for  a1=lowd1:uppd1
                CS_Mask(ia1,ia2,ia3)=Elliptic_Mask(ia1,ia2,ia3);
                if( (CorrFact/(CompSensingAcc*RelFill)*ProbMap(ia1,ia2,ia3)<(rand(1)) ) && Elliptic_Mask(ia1,ia2,ia3)==1)
                    CS_Mask(ia1,ia2,ia3)=0;
                    CSMask_Filling=CSMask_Filling-1;
                end
                
                ia1=ia1+1;
            end
            ia2=ia2+1;
        end
        ia3=ia3+1;
    end
    
    
    if (GoalFill > CSMask_Filling)CorrFact=CorrFact*1.05;
    else CorrFact=CorrFact/1.03;
    end
    
%      stepCS
%     CorrFact
%      GoalFill
%      CSMask_Filling
%abs(GoalFill-CSMask_Filling)*1.0/GoalFill
end

% }//if(m_bIsCompSensing){

US_MASK=fftshift(fftshift(fftshift(CS_Mask,1),2),3);
%US_MASK=fftshift(fftshift(US_MASK,1),2);
%imagesc(kmask)
%US_MASK=US_MASK;
% hist(Mask_DistRel(CS_Mask>0))




%{
        kmask=fftshift(fftshift(fftshift(kmask,1),2),3);
        Size_grid=size(kmask);
        
        for a=1:Size_grid(1);
            for b=1:Size_grid(2);
                for c=1:Size_grid(3);
                    % radius(a,b,c)=sqrt((a-Size_grid(1)*0.5-1)^2+ (b-Size_grid(2)*0.5-1)^2 +(c-Size_grid(3)*0.5-1)^2)+1;
                    radius(a,b,c)=(a-Size_grid(1)*0.5-1)^2 + (b-Size_grid(2)*0.5-1)^2 +(c-Size_grid(3)*0.5-1)^2 -safe_radius^2;
                    if radius(a,b,c)<0;radius(a,b,c)=0;else
                        radius(a,b,c)=sqrt( radius(a,b,c)/((Size_grid(1)*0.5*Size_grid(2)*0.5*Size_grid(3)*0.5)^(2/3)-safe_radius^2) );
                    end
                end
            end
        end
        
        %Prob_Map=(radius.^(expo));
        
        Prob_Map=(1-radius).^expo;
        
        if expo==0
            Prob_Map=double(radius==0);
        end
        
        %US_MASK=fftshift(fftshift(kmask,1),2);
        US_MASK=kmask;
        Tot_AcqP=sum(kmask(:));
        fil=1;
        while(fil>Fil_Fact)
            a=round(rand()*(Size_grid(1)-1))+1;
            b=round(rand()*(Size_grid(2)-1))+1;
            c=round(rand()*(Size_grid(3)-1))+1;
            if(Prob_Map(a,b,c)<rand()) ;
                US_MASK(a,b,c)=0;
            end
            fil=sum(US_MASK(:))/Tot_AcqP;
        end
        US_MASK=fftshift(fftshift(fftshift(US_MASK,1),2),3);
        %US_MASK=fftshift(fftshift(US_MASK,1),2);
        %imagesc(kmask)
        %figure
        %imagesc(US_MASK)
    end
%}
