function [total, injected, noise]=generatefakedatalam(lams,parangs, imdim,rpois,radii,sigstr)

nplanets=length(radii);

mu = [imdim*1.5 imdim*1.5];
Sigma = [5*imdim 0; 0 5*imdim];
x1 = 1:1:3*imdim;
x2 = 1:1:3*imdim;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));

base=F;

for i=1:3*imdim
    for j=1:3*imdim
        if sqrt((i-ceil(3*imdim/2))^2+(j-ceil(3*imdim/2))^2)<imdim/14
            base(i,j)=0;
        end
    end
end

base=base/max(max(base));
parangs=parangs*pi/180;

reflam=lams(ceil(length(lams)/2));

for i=1:length(parangs)
    for j=1:length(lams) 
        pos = length(lams)*(i-1)+j;
        scaledimage = imresize(base,lams(j)/reflam);
        [scalesize, dummy]= size(scaledimage);
        if not(rem(scalesize,2)) & scalesize > 3*imdim %%even and bigger
            tempres = imtranslate(scaledimage,[-.5 -.5]);
            croppedimages(:,:,pos) = imcrop(tempres,[1 2 3],[ceil(scalesize-3*imdim)/2 ceil(scalesize-3*imdim)/2 3*imdim-1  3*imdim-1]);
        elseif scalesize > 3*imdim %%odd and bigger
            croppedimages(:,:,pos) = imcrop(scaledimage,[1 2 3],[(scalesize-3*imdim)/2+1 (scalesize-3*imdim)/2+1 3*imdim-1  3*imdim-1]);
        elseif not(rem(scalesize,2)) %%even and smaller
            tempres = imtranslate(scaledimage,[.5 .5],'OutputView','full');
            croppedimages(:,:,pos) = zeropad(tempres,floor((3*imdim-scalesize)/2));
        else 
            croppedimages(:,:,pos) = zeropad(scaledimage,floor((3*imdim-scalesize)/2));
        end
    end
end
total=croppedimages(imdim+1:2*imdim,imdim+1:2*imdim,:);


thetas=[pi/3 2*pi/3 0];

for k=1:length(parangs)
    for n=1:length(lams)
        pos = length(lams)*(k-1)+n;
        for m=1:nplanets
            if m<4
                i(m)=ceil(imdim/2)+radii(m)*cos(thetas(m)+parangs(k));
                j(m)=ceil(imdim/2)+radii(m)*sin(thetas(m)+parangs(k));
            else
                i(m)=ceil(imdim/2)-radii(m)*cos(thetas(m-3)+parangs(k));
                j(m)=ceil(imdim/2)-radii(m)*sin(thetas(m-3)+parangs(k));
            end
        end

        [inject]=injectplanets(total(:,:,pos),i,j,imdim,sigstr);
        injected(:,:,pos)=inject-total(:,:,pos);
        specked(:,:,pos)=addspeckle(inject,imdim,1-rpois);
        pois(:,:,pos)=imnoise(specked(:,:,pos),'poisson');
        total(:,:,pos)=imnoise(pois(:,:,pos),'speckle',.02*rpois);
    end
end

noise = total - injected;


end