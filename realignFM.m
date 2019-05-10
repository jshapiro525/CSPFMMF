function [stack] = realignFM(images,parangs,lams)


N = length(parangs);
reflam = lams(ceil(length(lams)/2));
[imdim, dummy] = size(images(:,:,1));

[rep_parangs, rep_lams] = repparlam(parangs,lams);

for i = 1:N

%     scaledimage = imresize(images(:,:,i),rep_lams(i)/reflam);
%     [scalesize, dummy]= size(scaledimage);
%     if not(rem(scalesize,2)) & scalesize > imdim %%even and bigger
%         tempres = imtranslate(scaledimage,[-.5 -.5]);
%         croppedimages(:,:,i) = imcrop(tempres,[1 2 3],[ceil(scalesize-imdim)/2 ceil(scalesize-imdim)/2 imdim-1  imdim-1]);
%     elseif scalesize > imdim %%odd and bigger
%         croppedimages(:,:,i) = imcrop(scaledimage,[1 2 3],[(scalesize-imdim)/2+1 (scalesize-imdim)/2+1 imdim-1  imdim-1]);
%     elseif not(rem(scalesize,2)) %%even and smaller
%         tempres = imtranslate(scaledimage,[.5 .5],'OutputView','full');
%         croppedimages(:,:,i) = zeropad(tempres,floor((imdim-scalesize)/2));
%     else 
%         croppedimages(:,:,i) = zeropad(scaledimage,floor((imdim-scalesize)/2));
%     end
%    aligned(:,:,i)=imrotate(croppedimages(:,:,i),rep_parangs(i),'bicubic','crop');
    aligned(:,:,i)=imrotate(images(:,:,i),rep_parangs(i),'bicubic','crop');

end
    
stack = sum(images,3);

end