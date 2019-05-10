close all; 
clear all; 
clc;

lams = 1:.1:1.7;
parangs = [1:4:21];

imdim = 61;
epsilon = 0.03;

[totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon);

pix = imdim^2; %total number of pixels in an image
n = length(totalpert(1,1,:)); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to

[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strins for parangs and lambdas (one for each image)

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
P = inv(sqrt(diag(lam,0)))*V';
Sh_un = P*(Css1)*P';

[Uh2,eigh2] = eig(Sh_un);
[eighvals_un,indsh] = sort(diag(real(eigh2)),1,'descend');
Uh2 = Uh2(:, indsh);


%% Computing FM
lammat = diag(lam,0);

cbarS = inv(sqrt(lammat))*V'*Css1*V*inv(sqrt(lammat));

% myfig(Sh_un)
% title('Actual, unperturbed covariance matrix')
% 
% myfig(cbarS)
% title('Estimated unperturbed covariance matrix')

myfig(Sh_un-cbarS)
title('Difference between actual and estimated unperturbed covariance matrix')

%% Calculatin rho, phi, psi, and nu

[nu,rho] = eig(cbarS);

[rho,inds] = sort(diag(rho),1,'descend');

nu = nu(:, inds);
rhomat = diag(rho);

figure()
plot((eigvals_un-rho)./eigvals_un*100)
title('Percent unperturbed Eigenvalue difference')

figure()
plot(eigvals_un)
title('calculated eigenvalues')

figure()
plot(rho)
title('estimated eigenvalues')

myfig(nu)
title('Estimated unperturbed Eigenvectors')

myfig(Uh2)
title('Calculated unperturbed Eigenvectors')
 
myfig((Uh2-nu)./Uh2*100)
title('Percent unperturbed Eigenvector difference')
% 
% pdifeigvectun = reshape(abs((Uh2-nu)./Uh2)*100,n^2,1);
% 
% bins = linspace(0, ceil(max(pdifeigvectun)), ceil(max(pdifeigvectun)));
% xc = histc(pdifeigvectun,bins);
% 
% figure()
% semilogx(bins,xc,'r')
% hold on
% h = line([median(pdifeigvectun) median(pdifeigvectun)],[0 xc(2)]);
% title('Distribution of percent error, unperturbed')
% ylabel('Number of entries in the eigenvectors')
% xlabel('Percent Error')
% legend('Distribution of Percent Error',horzcat(['Median = ' num2str(median(pdifeigvectun)) '%']))


