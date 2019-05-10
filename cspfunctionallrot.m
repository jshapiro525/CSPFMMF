function [CSP, eigvals, Z, K, croppedimages,total2] = cspfunctionallrot(Total, parangs, lams, imdim, modes)

p=imdim^2; %total number of pixels in an image
n = length(Total(1,1,:)); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to

[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strins for parangs and lambdas (one for each image)

for i=1:n  % scaling and rotating for each image
        
    scaledimage = imresize(Total(:,:,i),reflam/rep_lams(i)); 
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

    total2(:,:,i)=imrotate(Total(:,:,i),rep_parangs(i),'bicubic','crop');
    X2(:,i)=reshape(total2(:,:,i),p,1);
    X1(:,i)=reshape(croppedimages(:,:,i),p,1);

end

X2=X2-mean(X2);
X1=X1-mean(X1);

X1=X1';
X2=X2';

Rh = X1*X1'*1/trace(X1*X1');
Rf = X2*X2'*1/trace(X2*X2');

R = Rh+Rf;

[U0,Sig] = eig(R);
[Sig,inds] = sort(diag(Sig),1,'descend');
Sig = diag(Sig,0);
U0 = U0(:, inds);

P = inv(Sig)^(.5)*U0';
P = real(P);
Sh = P*Rh*P';
Sf = P*Rf*P';

[Uh,eigh] = eig(Sh);
[eighvals,indsh] = sort(diag(real(eigh)),1,'descend');
Uh = real(Uh(:, indsh));

eigvals = eighvals;

W = Uh'*P;
Z = W*X1;
G = sqrt(inv(diag(eighvals)));  %% Normalization
Z = G*Z;

K = n+1-modes;

ZK=Z(K:end,:);
temp=(ZK'*ZK);

for i=1:n
    sig=temp*X1(i,:)';
    CSP(:,:,i)=reshape(sig',imdim,imdim);
end

% ZK = Z(1:K-1,:);
% temp=(ZK'*ZK);
% 
% for i=1:n
%     sig=(eye(imdim^2)-temp)*X1(i,:)';
%     CSP(:,:,i)=reshape(sig',imdim,imdim);
% end


end
