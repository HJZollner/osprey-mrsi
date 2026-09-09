function SparseData_ckkk = FillKDataWithRandom(SparseData_ckkk,kmask_kkk )
KSize=size(kmask_kkk);
[X,Y,Z] = ndgrid(1:KSize(1), 1:KSize(2),1:KSize(3));
xc=floor(KSize(1)/2)+1;
yc=floor(KSize(2)/2)+1;
zc=floor(KSize(3)/2)+1;
temp = ((X-xc).^2 + (Y-yc).^2+ (Z-zc).^2).^0.5 ;
Radius = 1+fftshift(fftshift(fftshift(temp,1),2),3);
for c=1:size(SparseData_ckkk,1)
    p = polyfit(log(Radius(kmask_kkk>0)),log(abs(SparseData_ckkk(c,kmask_kkk>0)))',1);
    AmpRandom=exp(p(1)*log(Radius)+p(2));
    RandData_kkk=AmpRandom.*(randn(KSize)+1j*randn(KSize));
    SparseData_ckkk(c,kmask_kkk==0)=RandData_kkk(kmask_kkk==0);
end
end

