close all; 
clear all; 
clc;

 lams = 1:.1:1.9;
% lams = 1;
parangs = [1:5:21];
% parangs = 1;
imdim = 71;
pix= imdim^2;
epsilon = 0.06;

n = length(parangs)*length(lams); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to
[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strins for parangs and lambdas (one for each image)
r = [18.5 18.5 18.5 14.5 20.5 14.5];
%
[totalpert2, injected2, noise2]=generatefakedatalam(lams,parangs,imdim,.05,r,epsilon);
totalpert = totalpert2(floor(imdim/2):end,floor(imdim/2):end,:);
injected = injected2(floor(imdim/2):end,floor(imdim/2):end,:);
noise = noise2(floor(imdim/2):end,floor(imdim/2):end,:);

imdim2 = length(noise(:,1,1));
pix2= imdim2^2;
%

[S12, S2, Css1, V, lammat, rhomat, nu, Z_un] = getunpert(noise2, reflam, rep_lams, rep_parangs, imdim);
[A12, dellammat, delrhomat, delV, delnu] = getpert(injected2, reflam, rep_lams, rep_parangs, imdim, lammat, rhomat, V, nu, S12, S2, Css1);


[S1, S2, Css1, V, lammat, rhomat, nu, Z_un] = getunpert(noise, reflam, rep_lams, rep_parangs, imdim2);
[A1, dellammat, delrhomat, delV, delnu] = getpert(injected, reflam, rep_lams, rep_parangs, imdim2, lammat, rhomat, V, nu, S1, S2, Css1);

delZ = rhomat^(-1/2) * nu' * lammat^(-1/2) * V' * A1 +...
       rhomat^(-1/2) * nu' * lammat^(-1/2) * delV' * S1 +...
       rhomat^(-1/2) * nu' * -1/2* lammat^(-3/2)*dellammat * V' * S1 +...
       rhomat^(-1/2) * delnu' * lammat^(-1/2) * V' * S1 + ...
       -1/2*rhomat^(-3/2)*delrhomat * nu' * lammat^(-1/2) * V' * S1;

[CSPsmall, eigvals, Z, K, croppedimages,total2] = cspfunctionallrot(totalpert, parangs, lams, imdim2, .01);


CSP = zeros(imdim,imdim,n);
CSP(floor(imdim/2):end,floor(imdim/2):end,:) = CSPsmall;

finalstack = realign(CSP,parangs,lams);
myfig(finalstack)
title('Result')

figure()
plot(eigvals)

X1 = A12+S12;

delZK = delZ(K:end,:);

for i = 1:length(delZK(:,1))
    pmodes(:,:,i) = reshape(delZK(i,:),imdim2, imdim2);
end

PMODES = zeros(imdim,imdim,n-K+1);
PMODES(floor(imdim/2):end,floor(imdim/2):end,:) = pmodes;

for i = 1:length(delZK(:,1))
    delzK(i,:) = reshape(PMODES(:,:,i),1,pix);
end

temp = (delzK'*delzK);
 
for i=1:n
    sig=temp*X1(i,:)';
    planetsig(:,:,i)=reshape(sig,imdim,imdim);
end


% delZK = delZ(1:K-1,:);
% temp=(delZK'*delZK);
% 
 
% for i=1:n
%     sig=(eye(pix)-temp)*X1(i,:)';
%     planetsig(:,:,i)=reshape(sig',imdim,imdim);
% end

planetstack = realign(planetsig,parangs,lams);

myfig(injected(:,:,1))
title('injected location')

myfig(planetstack)
title('Forward Model of planet signal')

% fitswrite(modes,'modes.fits')
% fitswrite(pmodes,'pmodes.fits')

