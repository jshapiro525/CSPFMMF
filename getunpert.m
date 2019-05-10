function [S1, S2, Css1, V, lammat, rhomat, nu, Z_un] = getunpert(noise, reflam, rep_lams, rep_parangs, imdim)

pix = imdim^2;
n = length(rep_lams);

for i=1:n  % scaling and rotating for each image
        
    scaledimage = imresize(noise(:,:,i),reflam/rep_lams(i)); 
    [scalesize, dummy]= size(scaledimage);
    
    % even and bigger:
    if not(rem(scalesize,2)) & scalesize > imdim 
        tempres = imtranslate(scaledimage,[-.5 -.5]);
        croppedimages(:,:,i) = imcrop(tempres,[1 2 3],[ceil(scalesize-imdim)/2 ceil(scalesize-imdim)/2 imdim-1  imdim-1]);
        
    % Odd and bigger:    
    elseif scalesize > imdim 
        croppedimages(:,:,i) = imcrop(scaledimage,[1 2 3],[(scalesize-imdim)/2+1 (scalesize-imdim)/2+1 imdim-1  imdim-1]);
        
    % Even and smaller:    
    elseif not(rem(scalesize,2))
        tempres = imtranslate(scaledimage,[.5 .5],'OutputView','full');
        croppedimages(:,:,i) = zeropad(tempres,floor((imdim-scalesize)/2));
        
    % Odd and Smaller:    
    else 
        croppedimages(:,:,i) = zeropad(scaledimage,floor((imdim-scalesize)/2));
    end

    total2(:,:,i)=imrotate(noise(:,:,i),rep_parangs(i),'bicubic','crop');
    S2(:,i)=reshape(total2(:,:,i),pix,1);
    S1(:,i)=reshape(croppedimages(:,:,i),pix,1);
    
end

S2=S2-mean(S2);
S1=S1-mean(S1);

S1=S1';
S2=S2';

Css1 = S1*S1'/(pix-1);%*1/trace(S1*S1');
Css2 = S2*S2'/(pix-1);%*1/trace(S2*S2');

[V,lam] = eig(Css1+Css2);
[lam,inds] = sort(diag(lam),1,'descend');
V = V(:, inds);
for k=1:n
    if V(1,k) < 0
        V(:,k) = V(:,k)*-1;
    end
end

P = inv(sqrt(diag(lam,0)))*V';
Sh_un = P*(Css1)*P';

[nu,rho] = eig(Sh_un);
[rho,indsh] = sort(diag(real(rho)),1,'descend');
nu = real(nu(:, indsh));
for k=1:n
    if nu(1,k) < 0
        nu(:,k) = nu(:,k)*-1;
    end
end

W_un = nu'*P;
Z_un = sqrt(inv(diag(rho)))*W_un*S1;

rhomat=diag(rho);
lammat=diag(lam);

end