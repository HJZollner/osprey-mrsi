function FreqMusic=DetermineSingleLowFreq(data_trrr,fMax,mrsiWaterParams)

%% Parameters
Fe = mrsiWaterParams.Water_mrProt.samplerate; % [Hz]
p = 3; % Number of pole pairs
N = size(data_trrr,1);

t = (0 : N-1) / Fe;


f = linspace (-fMax, fMax, 1024);
origf = linspace (-Fe/2, Fe/2, N);

FreqMusic=zeros(size(data_trrr,2),size(data_trrr,3),size(data_trrr,4));


for a=1:size(data_trrr,2)
    for b=1:size(data_trrr,3)
        parfor c=1:size(data_trrr,4)
        xn=data_trrr(:,a,b,c);
        %% Music

        [S, freq] = pmusic (xn, p, f, Fe);
        [Max, IxMax] = max (S);
      %  if( IxMax==numel(S)) | (IxMax==0) % estimaion is inaccurate
        %	FreqMusic(a,b,c) =0;
       % else
        	FreqMusic(a,b,c) = freq(IxMax);
	%end
    end
end


end
