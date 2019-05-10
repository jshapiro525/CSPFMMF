function [ n, medianerror, meanerror] = testFMimagesshort(lams,parangs)

epsilon = 0.03;
modes = 2;


imdim = 71;
[totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon);

%% Getting actual cov matrix and Z from perturbed data

pix = imdim^2; %total number of pixels in an image
n = length(totalpert(1,1,:)); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to

for i=1:n
    totalpert(:,:,i) = totalpert(:,:,i);
    noise(:,:,i) = noise(:,:,i);
end

[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strins for parangs and lambdas (one for each image)

X1 = [];
X2 = [];

for i=1:n  % scaling and rotating for each image
        
    scaledimage = imresize(totalpert(:,:,i),reflam/rep_lams(i)); 
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

    total2(:,:,i)=imrotate(totalpert(:,:,i),rep_parangs(i),'bicubic','crop');
    X2(:,i)=reshape(total2(:,:,i),pix,1);
    X1(:,i)=reshape(croppedimages(:,:,i),pix,1);

end

X2=X2-mean(X2);
X1=X1-mean(X1);

X1=X1';
X2=X2';

C1 = X1*X1'/(pix-1);
C2 = X2*X2'/(pix-1);;%*1/trace(X2*X2');

C = C1+C2;

[U0,Sig] = eig(C);
[Sig,inds] = sort(diag(Sig),1,'descend');
Sig = diag(Sig,0);
U0 = U0(:, inds);

for k=1:n
    if U0(1,k) < 0
        U0(:,k) = U0(:,k)*-1;
    end
end

P = inv(sqrt(Sig))*U0';
P = real(P);
Sh = P*C1*P';

[Uh1,eigh1] = eig(Sh);
[eighvals1,indsh] = sort(diag(real(eigh1)),1,'descend');
Uh1 = real(Uh1(:, indsh));

for k=1:n
    if Uh1(1,k) < 0
        Uh1(:,k) = Uh1(:,k)*-1;
    end
end

eigvals_true = eighvals1;

W_true = Uh1'*P;
Z_true = sqrt(inv(diag(eigvals_true)))*W_true*X1;

%% unperturbed Cov, eigenvals, and Z

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

[Uh2,eigh2] = eig(Sh_un);
[eighvals2,indsh] = sort(diag(real(eigh2)),1,'descend');
Uh2 = real(Uh2(:, indsh));
for k=1:n
    if Uh2(1,k) < 0
        Uh2(:,k) = Uh2(:,k)*-1;
    end
end

eigvals_un = eighvals2;

W_un = Uh2'*P;
Z_un = sqrt(inv(diag(eigvals_un)))*W_un*S1;

%% Signals and their perturbations

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

%% Building P
for k=1:n
    dellam(k)=V(:,k)'*Cas*V(:,k);
end

Gamma = lam + dellam';

for k=1:n
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

U = V + delV;

P = inv(sqrt(diag(Gamma,0)))*U';

%% % Check U

eigvect1error = abs((U0-U)./U0)*100;


%% Computing FM
Gammamat = diag(Gamma,0);
lammat = diag(lam,0);
dellammat = diag(dellam,0);

cbarS = inv(sqrt(lammat))*V'*Css1*V*inv(sqrt(lammat));

a = -3/2;
cbarA = inv(sqrt(Gammamat))*(V'*Css1*delV + delV'*Css1*V + U'*Cas*U)*inv(sqrt(Gammamat)) - 1/2*inv(sqrt(lammat))*V'*Css1*V*lammat^a*dellammat - 1/2*lammat^a*dellammat*V'*Css1*V*inv(sqrt(lammat));
Cbar = cbarS+cbarA;

%% Checking Cbar

cbarerror = abs((Cbar-Sh)./Sh)*100; 


temp = reshape(cbarerror,n*n,1);

meanerror = mean(temp);
medianerror = median(temp);

end