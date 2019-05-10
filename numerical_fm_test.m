% close all; 
% clear all; 
% clc;
% 
epsilon = 0.01;
modes = 2;
% 
% lams = 1.2:.2:1.8;
% parangs = [1:5:21];
% imdim = 71;
% [totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon);
% 
% 
imdim = 75;
load('betpicandinjections75_1.mat')
injected = A*epsilon;
totalpert = noise + injected;

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
    k/n*100
end

U = V + delV;

P = inv(sqrt(diag(Gamma,0)))*U';
%% Intermediate Checking

%% % Check C1
% myfig(C1)
% title('Actual C1')
% myfig(Css1+CAS1)
% title('Estimated C1')
% myfig(C1-(Css1+CAS1))
% title('Raw C1 difference')
% myfig(abs(C1-(Css1+CAS1))./C1*100)
% title('Percent C1 Difference')
% graph_element_distro(abs(C1-(Css1+CAS1))./C1*100,'Percent C1 Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

 
%% % Check C2
% myfig(C2)
% title('Actual C2')
% myfig(Css2+CAS2)
% title('Estimated C2')
% myfig(C2-(Css2+CAS2))
% title('Raw C2 difference')
% myfig(abs(C2-(Css2+CAS2))./C2*100)
% title('Percent C2 Difference')

%% Check Gamma
% 
% figure()
% plot(1:n,diag(Sig),'b',1:n,Gamma,'r')
% legend('True Gamma','Estimated Gamma')
% title('Eigenvalue Comparison')
% 
% figure()
% plot(diag(Sig)-Gamma)
% title('Raw Eigenvalue Difference')
% 
figure()
plot(abs(diag(Sig)-Gamma)./diag(Sig)*100)
title('Percent Eigenvalue Difference')

%% % Check U
% 
% 
myfig(U0)
title('Actual U')
myfig(U)
title('Estimated U')
myfig(U0-U)
title('Raw Eigenvector differences')
% myfig(abs(U0-U)./U0*100)
% title('Percent Eigenvector Difference')

pdifeig= abs((U0-U)./U0)*100;

graph_element_distro(pdifeig,'Percent Eigenvector Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

%% Do Gamma and U work with C1?

% P_est = inv(sqrt(diag(Gamma,0)))*U';
% cbar_est  = P*C1*P';
% 
% myfig(Sh)
% title('Actual Cbar')
% 
% myfig(cbar_est)
% title('Estimated Cbar from U, gamma and C1')
% 
% myfig(cbar_est-Sh)
% title('Raw difference')
% 
% myfig((cbar_est-Sh)./Sh*100)
% title('percent difference')
% 
 graph_element_distro(abs((cbar_est-Sh)./Sh)*100,'Percent Cbar Difference','Percent Error','Number of entries in U','Distribution of Percent Error')


%% Computing FM
Gammamat = diag(Gamma,0);
lammat = diag(lam,0);
dellammat = diag(dellam,0);

cbarS = inv(sqrt(lammat))*V'*Css1*V*inv(sqrt(lammat));

a = -3/2;
cbarA = inv(sqrt(Gammamat))*(V'*Css1*delV + delV'*Css1*V + U'*Cas*U)*inv(sqrt(Gammamat)) - 1/2*inv(sqrt(lammat))*V'*Css1*V*lammat^a*dellammat - 1/2*lammat^a*dellammat*V'*Css1*V*inv(sqrt(lammat));
Cbar = cbarS+cbarA;

%% Checking Cbar
% % 
% % myfig(Sh)
% % title('Actual Cbar')
% 
% % myfig(Cbar)
% % title('Estimated, FM Cbar')
% 
% % myfig(abs(Cbar-Sh))
% % title('Full Cbar Raw difference')
% 
% perdifcbar = abs((Cbar-Sh)./Sh)*100; 
% pdifcbartoun = abs((cbarS-Sh)./Sh)*100;
% 
% myfig(perdifcbar)
% title('Full Cbar percent difference')
% 

graph_element_distro(perdifcbar,'Percent FM Cbar Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

%% Computing Z FM

%% Calculatin rho, phi, psi, and nu

[nu,rho] = eig(cbarS);

[rho,inds] = sort(diag(rho),1,'descend');

nu = nu(:, inds);
for k=1:n
    if nu(1,k) < 0
        nu(:,k) = nu(:,k)*-1;
    end
end

rhomat = diag(rho);

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
Phi = rho + delrho';
Psi = nu + delnu;

%% Intermediate Checking
% figure()
% plot(eigvals_true)
% hold on
% plot(Phi)
% legend('True Eigenvalues','Estimated Eigenvalues')
% 
figure()
plot(abs((eigvals_true-Phi)./eigvals_true*100))
title('Percent Eigenvalue difference')
% 
% figure()
% plot(abs((eigvals_un-rho)./eigvals_un*100))
% title('Percent unperturbed Eigenvalue difference')
% 
% figure()
% plot(abs((eigvals_true-rho)./eigvals_un*100))
% title('Percent Eigenvalue difference (true to unpert.)')
% 
% myfig(Uh1)
% title('True Eigenvectors')
% 
% myfig(Psi)
% title('Estimated Eigenvectors')
% 
% myfig(abs((Uh1-Psi)./Psi)*100)
% title('Percent Eigenvector difference')

myfig(abs((Uh2-nu)./Uh2)*100)
title('Percent unperturbed Eigenvector difference')

pdifeigvectun = abs((Uh1-nu)./Uh2)*100;
graph_element_distro(pdifeigvectun,'Distribution of perturbed eigenvector percent error to unperturbed model','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');


pdifeigvect = abs((Uh1-Psi)./Psi)*100;
graph_element_distro(pdifeigvect,'Distribution of Perturbed Eigenvector percent error','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');


%% Calculate Z

Zs = rhomat^(-1/2)*nu'*lammat^(-1/2)*V'*S1;
delZ = rhomat^(-1/2) * nu' * lammat^(-1/2) * V' * A1 +...
       rhomat^(-1/2) * nu' * lammat^(-1/2) * delV' * S1 +...
       rhomat^(-1/2) * nu' * -1/2* lammat^(-3/2)*dellammat * V' * S1 +...
       rhomat^(-1/2) * delnu' * lammat^(-1/2) * V' * S1 + ...
       -1/2*rhomat^(-3/2)*delrhomat * nu' * lammat^(-1/2) * V' * S1;
   
Z = Zs + delZ;

myfig(Z_un-Zs)
title('Difference in unperturbed projection matrices')

pdifZ = abs((Z_true-Z)./Z_true)*100;
pdifZun = abs((Z_true-Zs)./Z_true)*100;
graph_element_distro(pdifZun,'Difference in perturbed projection matrices to unperturbed model','Percent Error','Number of entries in Z','Distribution of Percent Error');

myfig(log10(pdifZ))
title('Log of Percent Difference in Z')

graph_element_distro(pdifZ,'Difference in perturbed projection matrices','Percent Error','Number of entries in Z','Distribution of Percent Error');

K=n-modes;

X1 = A1+S1;

delZK = delZ(K:end,:);
temp=(delZK'*delZK);
temp2 = delZ'*delZ;
 
for i=1:n
    sig=temp*X1(i,:)';
    sig2=temp2*X1(i,:)';
    planetsig(:,:,i)=reshape(sig,imdim,imdim);
    fullpsig(:,:,i)=reshape(sig2,imdim,imdim);
end

planetstack = realign(planetsig,parangs,lams);
fullpstack = realign(fullpsig,parangs,lams);

myfig(injected(:,:,1))
title('injected location')

myfig(planetstack)
title('Forward Model of planet signal')

myfig(fullpstack)
title('Full Forward Model of planet signal')

% Ztemp = inv(sqrt(rhomat))*Ws*S1;


% myfig(Z_un)
% title('Actual, unperturbed Z matrix')
% % 
% myfig(Zs)
% title('Estimated unperturbed Z matrix')
% 
% myfig(Z_un-Zs)
% title('Difference between actual and estimated unperturbed Z matrix')
% 
% myfig(Z_true)
% title('Actual Z Matrix')
% 
% myfig(Z)
% title('Estimated FM Z')
% 
% myfig(abs(Z_true-Z))
% title('Z Raw difference')
% 
% perdifZ = abs((Z_true-Z)./Z_true)*100; 
% pdifZtoun = abs((Z_true-Zs)./Z_true)*100;
% 
% myfig(perdifZ)
% title('Full Z percent difference')
% 
% perdifZ = reshape(perdifZ,n*pix,1);
% pdifZtoun = reshape(pdifZtoun,n*pix,1);
% 
% bins = linspace(0, ceil(max(perdifZ))/10000, ceil(max(perdifZ)));
% xc = histc(perdifZ,bins);
% 
% figure()
% semilogx(bins,xc,'r')
% hold on
% h = line([median(perdifZ) median(perdifZ)],[0 xc(2)]);
% title('Distribution of percent error')
% ylabel('Number of entries in Z')
% xlabel('Percent Error')
% legend('Distribution of Percent Error',horzcat(['Median = ' num2str(median(perdifcbar)) '%']))
% 
% bins = linspace(0, ceil(max(pdifZtoun)), ceil(max(pdifZtoun)));
% xc = histc(pdifZtoun,bins);
% 
% figure()
% semilogx(bins,xc,'r')
% hold on
% h = line([median(pdifZtoun) median(pdifZtoun)],[0 xc(2)]);
% title('Distribution of percent error to just the unperturbed model')
% ylabel('Number of entries in Z')
% xlabel('Percent Error')
% legend('Distribution of Percent Error',horzcat(['Median = ' num2str(median(pdifcbartoun)) '%']))
