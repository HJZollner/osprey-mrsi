function h = volimagesc( Volume,sc )
if nargin==2
    h = imagesc(Vol2Image(Volume),sc);
elseif nargin==1
     h = imagesc(Vol2Image(Volume));
end

end

