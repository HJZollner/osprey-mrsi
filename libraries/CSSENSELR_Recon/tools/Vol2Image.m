function [plotImage] = Vol2Image( Vol )

plotImage= zeros(1, size(Vol,2)*ceil(sqrt(size(Vol,3))));
z=1;
while z<=size(Vol,3)
    TempImage= Vol(:,:,z);
    z=z+1;
    for ind=2:ceil(sqrt(size(Vol,3)))
        if z<=size(Vol,3)
            TempImage =[ TempImage Vol(:,:,z) ];
        else
            TempImage =[ TempImage 0*Vol(:,:,1) ];
        end
        z=z+1;
    end
    plotImage=[plotImage;TempImage];
end

plotImage=plotImage(2:end,:);
end

