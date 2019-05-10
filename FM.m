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
r = [20.5 20.5 20.5 26.5 12.5 32.5];

n = length(parangs)*length(lams); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to
[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strins for parangs and lambdas (one for each image)

[totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,r,epsilon);

[S1, S2, Css1, V, lammat, rhomat, nu, Z_un] = getunpert(noise, reflam, rep_lams, rep_parangs, imdim);
[A1, dellammat, delrhomat, delV, delnu] = getpert(injected, reflam, rep_lams, rep_parangs, imdim, lammat, rhomat, V, nu, S1, S2, Css1);

Zun = rhomat^(-1/2) * nu' * lammat^(-1/2) * V' * S1;

delZ = rhomat^(-1/2) * nu' * lammat^(-1/2) * V' * A1 +...
       rhomat^(-1/2) * nu' * lammat^(-1/2) * delV' * S1 +...
       rhomat^(-1/2) * nu' * -1/2* lammat^(-3/2)*dellammat * V' * S1 +...
       rhomat^(-1/2) * delnu' * lammat^(-1/2) * V' * S1 + ...
       -1/2*rhomat^(-3/2)*delrhomat * nu' * lammat^(-1/2) * V' * S1;

[CSP, eigvals, Z, K, croppedimages,total2] = cspfunctionallrot(totalpert, parangs, lams, imdim, .01);

finalstack = realign(CSP,parangs,lams);
myfig(finalstack)
title('Result')

figure()
plot(eigvals)

X1 = A1+S1;

Zest = delZ + Zun;

delZK = delZ(K:end,:);
ZestK = Zest(K:end,:);

temp = (delZK'*delZK);
temp2 = delZ'*delZ;
temp3 = ZestK'*ZestK; 
 
for i=1:n
    sig=temp*X1(i,:)';
    sig2=temp2*X1(i,:)';
    sig3 = temp3*X1(i,:)';
    planetsig(:,:,i)=reshape(sig,imdim,imdim);
    fullpsig(:,:,i)=reshape(sig2,imdim,imdim);
    fullsigest(:,:,i)=reshape(sig3,imdim,imdim);
end


% delZK = delZ(1:K-1,:);
% temp=(delZK'*delZK);
% 
 
% for i=1:n
%     sig=(eye(pix)-temp)*X1(i,:)';
%     planetsig(:,:,i)=reshape(sig',imdim,imdim);
% end

planetstack = realign(planetsig,parangs,lams);
fullpstack = realign(fullpsig,parangs,lams);
fulleststack = realign(fullsigest,parangs,lams);

myfig(injected(:,:,1))
title('injected location')

myfig(planetstack)
title('Forward Model of planet signal')

poscomparison = planetstack+injected(:,:,1)*10^4;
myfig(poscomparison)
title('Transformations compared to injected')
colormap('hsv')

myfig(fulleststack)
title('Complete Forward Model')

myfig(fullpstack)
title('Full Forward Model of planet signal')

for i = 1:n
    modes(:,:,i) = reshape(Z(i,:),imdim, imdim);
    pertmodes(:,:,i) = reshape(delZ(i,:),imdim, imdim);
end

% a = finalstack(41:51,41:51);
% b = fulleststack(41:51,41:51);
% a = a/(mean(mean(a)));
% b = b/(mean(mean(b)))*max(max(a))/(max(max(b)));
% 
% figure()
% plot([1:11],a(7,:),[1:11],b(7,:))
% legend('Science Value',' FM Estimate')
% xlabel('pixel')
% ylabel('value')
% 
% figure()
% plot([1:11],a(:,6),[1:11],b(:,6))
% legend('Science Value',' FM Estimate')
% xlabel('pixel')
% ylabel('value')

