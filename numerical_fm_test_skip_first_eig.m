close all; 
clear all; 
clc;

epsilon = 0.03;
modes = 18;

% lams = 1.2:(1.8-1.2)/36:1.8;
% parangs = [1:31/36:32];
% imdim = 71;
% [totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon);


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

P_true = inv(sqrt(Sig))*U0';
P_true = real(P_true);
Sh = P_true*C1*P_true';

[Uh1,eigh1] = eig(Sh);
[eighvals1,indsh] = sort(diag(real(eigh1)),1,'descend');
Uh1 = real(Uh1(:, indsh));

for k=1:n
    if Uh1(1,k) < 0
        Uh1(:,k) = Uh1(:,k)*-1;
    end
end

eigvals_true = eighvals1;

W_true = Uh1'*P_true;
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

P_un = inv(sqrt(diag(lam,0)))*V';
Sh_un = P_un*(Css1)*P_un';

[Uh2,eigh2] = eig(Sh_un);
[eighvals2,indsh] = sort(diag(real(eigh2)),1,'descend');
Uh2 = real(Uh2(:, indsh));
for k=1:n
    if Uh2(1,k) < 0
        Uh2(:,k) = Uh2(:,k)*-1;
    end
end

eigvals_un = eighvals2;

W_un = Uh2'*P_un;
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

Cest = Css1 + Css2 + Cas + (A1*A1'+A2*A2')/(pix-1);

[U, Gamma] = eig(Cest);
[Gamma,indest] = sort(diag(real(Gamma)),1,'descend');
U = real(U(:, indest));
for k=1:n
    if U(1,k) < 0
        U(:,k) = U(:,k)*-1;
    end
end

dellam = Gamma - lam;

delV = U-V;

P_est = inv(sqrt(diag(Gamma,0)))*U';
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
% log_graph_element_distro(abs(C1-(Css1+CAS1))./C1*100,'Percent C1 Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

 
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
% figure()
% plot(abs(diag(Sig)-Gamma)./diag(Sig)*100)
% title('Percent Eigenvalue Difference')

%% % Check U
% 
% 
% myfig(U0)
% title('Actual U')
% myfig(U)
% title('Estimated U')
% myfig(U0-U)
% title('Raw Eigenvector differences')
% % myfig(abs(U0-U)./U0*100)
% % title('Percent Eigenvector Difference')
% 
% pdifeig= abs((U0-U)./U0)*100;
% 
% log_graph_element_distro(pdifeig,'Percent Eigenvector Difference','Percent Error','Number of entries in U','Distribution of Percent Error')


%% Computing Cbar
Gammamat = diag(Gamma,0);
lammat = diag(lam,0);
dellammat = diag(dellam,0);

cbarS = inv(sqrt(lammat))*V'*Css1*V*inv(sqrt(lammat));

cbarA = P_est*(Css1+CAS1)*P_est'-cbarS;
Cbar = cbarS+cbarA;

%% Checking Cbar
% % 
myfig(Sh)
title('Actual Cbar')

myfig(Cbar)
title('Estimated, FM Cbar')

myfig(abs(Cbar-Sh))
title('Full Cbar Raw difference')

perdifcbar = abs((Cbar-Sh)./Sh)*100; 
pdifcbartoun = abs((cbarS-Sh)./Sh)*100;

myfig(perdifcbar)
title('Full Cbar percent difference')

log_graph_element_distro(perdifcbar,' ','Percent Error','Number of entries in U','Distribution of Percent Error')

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
    k/n*100
end

delrhomat = diag(delrho);
Phi = rho + delrho';
Psi = nu + delnu;

%% Intermediate Checking
figure()
plot(eigvals_true)
hold on
plot(Phi)
legend('True Eigenvalues','Estimated Eigenvalues')
% 
figure()
plot(abs((eigvals_true-Phi)./eigvals_true*100))
title('Percent Eigenvalue difference')
% 
% figure()
% plot(abs((eigvals_un-rho)./eigvals_un*100))
% title('Percent unperturbed Eigenvalue difference')

figure()
plot(abs((eigvals_true-rho)./eigvals_un*100))
title('Percent Eigenvalue difference (true to unpert.)')

myfig(Uh1)
title('True Eigenvectors')

myfig(Psi)
title('Estimated Eigenvectors')

myfig(Uh1-Psi)

myfig(abs((Uh1-Psi)./Psi)*100)
title('Percent Eigenvector difference')

myfig(abs((Uh2-nu)./Uh2)*100)
title('Percent unperturbed Eigenvector difference')

pdifeigvectun = abs((Uh1-nu)./Uh2)*100;
log_graph_element_distro(pdifeigvectun,'Distribution of perturbed eigenvector percent error to unperturbed model','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');


pdifeigvect = abs((Uh1-Psi)./Psi)*100;
log_graph_element_distro(pdifeigvect,'Distribution of Perturbed Eigenvector percent error','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');


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
log_graph_element_distro(pdifZun,'Difference in perturbed projection matrices to unperturbed model','Percent Error','Number of entries in Z','Distribution of Percent Error');

myfig(log10(pdifZ))
title('Log of Percent Difference in Z')

log_graph_element_distro(pdifZ,' ','Percent Error','Number of entries in Z','Distribution of Percent Error');

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

%% comparison after here
close all;

for modes =[1 2:2:10 13:3:37 40:5:75]
    
    K = n- modes; 
    %ZK = Z(K:end,:);
    Z_trueK = Z_true(1:modes,:);
    %temp1 = ZK'*ZK;
    temp_true= Z_trueK'*Z_trueK;

    for i=1:n
     %   sig=temp1*X1(i,:)';
        sig_true=(eye(pix)-temp_true)*X1(i,:)';
     %   res(:,:,i)=reshape(sig,imdim,imdim);
        restrue(:,:,i)=reshape(sig_true,imdim,imdim);
      %  i/n*100
    end

    %finres = realign(res,parangs,lams);
    finres_true = realignscale(restrue,parangs,lams);

    %myfig(finres)
    %title('predicted FM')
    myfig(finres_true)
    title(horzcat(['actual FM ' num2str(modes) ' modes']))
end

