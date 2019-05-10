function [mf] = fmmf(template, temploc, base, stampsize, modes)

[n,pix,nlams] = size(template);
y0=temploc(1);
x0=temploc(2);

imdim= sqrt(pix);
stampoff = floor(stampsize/2);

xc = ceil(imdim/2);
yc = ceil(imdim/2);
r0 = sqrt((x0-xc)^2+(y0-yc)^2);
th0 = atan2(yc-y0,xc-x0);

nmodes=length(modes);

for lams = 1:nlams
    for i=1:nmodes
        k = modes(i);
        templatesquare(:,:,i,lams) = reshape(template(k,:,lams),imdim,imdim);
        templatestamp(:,:,i,lams) = templatesquare(y0-stampoff:y0+stampoff,x0-stampoff:x0+stampoff,i);
        basesquare(:,:,i,lams) = reshape(base(k,:,lams),imdim,imdim);
    end

    for k=1:nmodes
        mf(:,:,k,lams) = zeros(imdim,imdim);
        for i = stampoff+1:imdim-stampoff-1
            for j = stampoff+1:imdim-stampoff-1
                x2 = j;
                y2 = i;

                r2 = sqrt((x2-xc)^2+(y2-yc)^2);
                th2 = atan2(yc-y2,xc-x2);
                scalefact=r2/r0;
                if scalefact == 0
                    continue
                end
                thr = th2-th0;

                scaledimage = imresize(templatestamp(:,:,k,lams),scalefact);
                [scalesize, dummy]= size(scaledimage);
                if not(rem(scalesize,2)) & scalesize > stampsize %%even and bigger
                    tempres = imtranslate(scaledimage,[-.5 -.5]);
                    croppedimagek = imcrop(tempres,[1 2 3],[ceil(scalesize-stampsize)/2 ceil(scalesize-stampsize)/2 stampsize-1  stampsize-1]);
                elseif scalesize > stampsize %%odd and bigger
                    croppedimagek= imcrop(scaledimage,[1 2 3],[(scalesize-stampsize)/2+1 (scalesize-stampsize)/2+1 stampsize-1  stampsize-1]);
                elseif not(rem(scalesize,2)) %%even and smaller
                    tempres = imtranslate(scaledimage,[.5 .5],'OutputView','full');
                    croppedimagek = zeropad(tempres,floor((stampsize-scalesize)/2));
                else 
                    croppedimagek = zeropad(scaledimage,floor((stampsize-scalesize)/2));
                end
                stamprot = imrotate(croppedimagek,thr*180/pi,'bicubic','crop');
                basestamp = basesquare(i-stampoff:i+stampoff,j-stampoff:j+stampoff,k,lams);

                mf(i,j,k,lams) = sum(sum(stamprot.*basestamp));
            end
        end
    end
end
        

end