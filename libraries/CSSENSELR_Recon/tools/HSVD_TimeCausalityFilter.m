function filtmrsiData_trrr  = HSVD_TimeCausalityFilter(mrsiData_trrr,Samplerate,modelOrder)
    filtmrsiData_trrr = zeros(size(mrsiData_trrr),class(mrsiData_trrr));% t-r-r-r
    parfor a=1:size(mrsiData_trrr,2)
        filtmrsiData_trrr(:,a,:,:)=InsideLoopFunction(squeeze(mrsiData_trrr(:,a,:,:)),Samplerate,modelOrder);
    end
end

function filtmrsiData_trr = InsideLoopFunction(mrsiData_trr,Samplerate,modelOrder)
    filtmrsiData_trr = zeros(size(mrsiData_trr),class(mrsiData_trr));% t-r-r
    for b=1:size(mrsiData_trr,2)
            for c=1:size(mrsiData_trr,3)
                InvTimeSerie=flip(squeeze(mrsiData_trr(:,b,c)));
                if sum(abs(InvTimeSerie(:)))>0
                    [~,damps,basis,amps] = HSVD(InvTimeSerie,Samplerate, modelOrder);
                    indx   = find((damps <=0) );
                    TimeSerie=mrsiData_trr(:,b,c)-flip(sum(basis(:, indx) * diag(amps(indx), 0), 2));
                    filtmrsiData_trr(:,b,c) = mrsiData_trr(:,b,c)-flip(sum(basis(:, indx) * diag(amps(indx), 0), 2));
                end	
        end
    end
end
