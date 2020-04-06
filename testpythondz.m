close all
clear all
clc

%%loading in, checking As

load('C:/Users/jshap/Documents/pyCSP/testmodeldata.mat')
n=size(models,1);%total number of images
imdim=281;
modelimgs = zeros(n,imdim^2);
modelimgs(:,inds) = models;

for i=1:n
    modelimgsreshape(:,:,i)=reshape(modelimgs(i,:),imdim,imdim);
end

fitswrite(modelimgsreshape,'modelimagesdata.fits')

pix = size(A1,2);

%% Computing FM

Css1 = S1*S1'/(pix-1);%*1/trace(S1*S1');
Css2 = S2*S2'/(pix-1);%*1/trace(S1*S1');

[Vcheck,lamcheck] = eig(Css1+Css2);
[lamcheck,indices] = sort(diag(lamcheck),1,'descend');
Vcheck = Vcheck(:, indices);

lamdiff=(lamcheck-diag(Lam))./lamcheck*100;


P_un = inv(sqrt(diag(lamcheck,0)))*Vcheck';
Cbar_un = P_un*(Css1)*P_un';

[Ycheck,Omegacheck] = eig(Cbar_un);
[Omegacheck,indsh] = sort(diag(real(Omegacheck)),1,'descend');
Ycheck = Ycheck(:, indsh);

omegadiff = (Omegacheck-diag(Omega))./Omegacheck*100;


CAS1= (A1*S1'+S1*A1')/(pix-1);
CAS2= (A2*S2'+S2*A2')/(pix-1);

Cas = CAS1+CAS2;
Css1 = S1*S1'/(pix-1);

% This part was to see how it is set.

Vcheck = V;
lamcheck = diag(Lam);
Omegacheck = diag(Omega);
Ycheck = Y;


for k=1:n
    dellam(k)=Vcheck(:,k)'*Cas*Vcheck(:,k);
end

Gamma = lamcheck + dellam';

for k=1:n
    tempvectp = zeros(n,1);
    for p=1:n
        if p==k
            continue
        else
            tempvectp=tempvectp+Vcheck(:,p)'*Cas*Vcheck(:,k)/(lamcheck(k)-lamcheck(p))*Vcheck(:,p);
        end
    end
    delV(:,k) = tempvectp;
end

U = Vcheck + delV;

Gammamat = diag(Gamma,0);
lammat = diag(lamcheck,0);
dellammat = diag(dellam,0);

cbarA = inv(sqrt(Gammamat))*(Vcheck'*Css1*delV + delV'*Css1*Vcheck + U'*Cas*U)*inv(sqrt(Gammamat)) - 1/2*inv(sqrt(lammat))*Vcheck'*Css1*Vcheck*lammat^(-3/2)*dellammat - 1/2*lammat^(-3/2)*dellammat*Vcheck'*Css1*Vcheck*inv(sqrt(lammat));

for k=1:n
    delomega(k)=Ycheck(:,k)'*cbarA*Ycheck(:,k);
    tempvectp = zeros(n,1);
    for p=1:n
        if p==k
            continue
        else
            tempvectp=tempvectp+Ycheck(:,p)'*cbarA*Ycheck(:,k)/(Omegacheck(k)-Omegacheck(p))*Ycheck(:,p);
        end
    end
    dely(:,k) = tempvectp;
end

delomegamat = diag(delomega);
omegamat= diag(Omegacheck);

delZ = omegamat^(-1/2) * Ycheck' * lammat^(-1/2) * Vcheck' * A1 +...
       omegamat^(-1/2) * Ycheck' * lammat^(-1/2) * delV' * S1 +...
       omegamat^(-1/2) * Ycheck' * -1/2* lammat^(-3/2)*dellammat * Vcheck' * S1 +...
       omegamat^(-1/2) * dely' * lammat^(-1/2) * Vcheck' * S1 + ...
       -1/2*omegamat^(-3/2)*delomegamat * Ycheck' * lammat^(-1/2) * Vcheck' * S1;
   
dZdiff=(delZ-dZ)./delZ*100;


%% plotting differences

figure()
yyaxis left
plot(diag(Lam),'rx','LineWidth',3)
hold on
plot(lamcheck,'b+','LineWidth',3)
yyaxis right
plot(lamdiff,'k-.','LineWidth',3)
legend('Python','MATLAB','Percent Difference','Fontsize',20)
xlabel('Eigenvalue','Fontsize',20)
ylabel('\% Difference','Fontsize',20)
hold off
yyaxis left
ylabel('Value','Fontsize',20)
set(gca,'FontSize',20)
title('Lambda eigenvalues percent difference')


figure()
yyaxis left
plot(diag(Omega),'rx','LineWidth',3)
hold on
plot(Omegacheck,'b+','LineWidth',3)
yyaxis right
plot(omegadiff,'k-.','LineWidth',3)
legend('Python','MATLAB','Percent Difference','Fontsize',20)
xlabel('Eigenvalue','Fontsize',20)
ylabel('\% Difference','Fontsize',20)
hold off
yyaxis left
ylabel('Value','Fontsize',20)
set(gca,'FontSize',20)
title('Omega eigenvalues percent difference')



dZfig=nan(n,imdim^2);
dZfig(:,inds) = dZ;
delZfig=nan(n,imdim^2);
delZfig(:,inds)=delZ;
Zactual=nan(n,imdim^2);
Zactual(:,inds)=Z;
for i=1:n
    dZfigresh(:,:,i)=reshape(dZfig(i,:),imdim,imdim);
    delZfigresh(:,:,i)=reshape(delZfig(i,:),imdim,imdim);
    Zresh(:,:,i) = reshape(Zactual(i,:),imdim, imdim);
end


fitswrite(delZfigresh,'dZ_Matlab.fits')
fitswrite(dZfigresh,'dZ_Python.fits')
fitswrite(Zresh,'Z_Python.fits')

mode = 1;
myfig(delZfigresh(:,:,mode))
title('MATLAB delta_z')
myfig(dZfigresh(:,:,mode))
title('Python delta_z')
myfig(delZfigresh(:,:,mode)-dZfigresh(:,:,mode))
title('absolute differences in delta_z')
myfig(Zresh(:,:,mode))
title('CSP Z mode')
% myfig((delZfigresh(:,:,mode)-dZfigresh(:,:,mode))./delZfigresh(:,:,mode)*100)
% title('Percent differences in delta_z')

% doesitmatch=sum(sum(dZ.*delZ))/sum(sum(delZ.*delZ));






