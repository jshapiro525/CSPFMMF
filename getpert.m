function [A1, dellammat, delrhomat, delV, delnu] = getpert(injected, reflam, rep_lams, rep_parangs, imdim, lammat, rhomat, V, nu, S1, S2, Css1)

pix = imdim^2;
n = length(rep_lams);

lam = diag(lammat);
rho = diag(rhomat);

for i=1:n  % scaling and rotating for each image
        
    scaledimage = imresize(injected(:,:,i),reflam/rep_lams(i)); 
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

    total2(:,:,i)=imrotate(injected(:,:,i),rep_parangs(i),'bicubic','crop');
    A2(:,i)=reshape(total2(:,:,i),pix,1);
    A1(:,i)=reshape(croppedimages(:,:,i),pix,1);
    
end

A2=A2-mean(A2);
A1=A1-mean(A1);

A1=A1';
A2=A2';

CAS1= (A1*S1'+S1*A1')/(pix-1);
CAS2= (A2*S2'+S2*A2')/(pix-1);

Cas = CAS1+CAS2;

for k=1:n
    
    dellam(k)=V(:,k)'*Cas*V(:,k);

    tempvectp = zeros(n,1);
    for p=1:n
        if p==k
            continue
        else
            tempvectp=tempvectp+V(:,p)'*Cas*V(:,k)/(lam(k)-lam(p))*V(:,p);
        end
    end
    delV(:,k) = tempvectp;
end

dellammat = diag(dellam);

Gammamat = lammat+dellammat;
U = V + delV;

cbarA = inv(sqrt(Gammamat))*(V'*Css1*delV + delV'*Css1*V + U'*Cas*U)*inv(sqrt(Gammamat)) - 1/2*inv(sqrt(lammat))*V'*Css1*V*lammat^(-3/2)*dellammat - 1/2*lammat^(-3/2)*dellammat*V'*Css1*V*inv(sqrt(lammat));


for k=1:n
    delrho(k)=nu(:,k)'*cbarA*nu(:,k);
    tempvectp = zeros(n,1);
    for p=1:n
        if p==k
            continue
        else
            tempvectp=tempvectp+nu(:,p)'*cbarA*nu(:,k)/(rho(k)-rho(p))*nu(:,p);
        end
    end
    delnu(:,k) = tempvectp;
end

delrhomat = diag(delrho);


end