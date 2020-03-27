load('C:/Users/jshap/Documents/pyCSP/testmodeldata.mat')
n=size(models,1);
imdim=281;
modelimgs = zeros(n,imdim^2);
modelimgs(:,inds) = models;

for i=1:n
    modelimgsreshape(:,:,i)=reshape(modelimgs(i,:),imdim,imdim);
end

fitswrite(modelimgsreshape,'modelimagesdata.fits')

totalpert = noise + injected;
pix = size(A1,1);
n = size(A1,1); %total number of images

Css1 = S1*S1'/(pix-1);%*1/trace(S1*S1');
Css2 = S2*S2'/(pix-1);%*1/trace(S2*S2');

[V,lam] = eig(Css1+Css2);
[lam,inds] = sort(diag(lam),1,'descend');
V = V(:, inds);

for k=1:n
    tempind = find(abs(V(:,k))==max(abs(V(:,k))));
    if sign(V(tempind,k)) ~= sign(U_true(tempind,k))
        V(:,k) = V(:,k)*-1;
    end
end

P_un = inv(sqrt(diag(lam,0)))*V';
Cbar_un = P_un*(Css1)*P_un';

[W_un,Phi_un] = eig(Cbar_un);
[Phi_un,indsh] = sort(diag(real(Phi_un)),1,'descend');
W_un = W_un(:, indsh);

for k=1:n
    tempind = find(abs(W_un(:,k))==max(abs(W_un(:,k))));
    if sign(W_un(tempind,k)) ~= sign(W_true(tempind,k))
        W_un(:,k) = W_un(:,k)*-1;
    end
end

W_un = W_un'*P_un;
Z_un = sqrt(inv(diag(Phi_un)))*W_un*S1;


CAS1= (A1*S1'+S1*A1')/(pix-1);
CAS2= (A2*S2'+S2*A2')/(pix-1);

Cas = CAS1+CAS2;

%% Building P and Cbar

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

P_true = inv(sqrt(diag(Gamma,0)))*U';

Gammamat = diag(Gamma,0);
lammat = diag(lam,0);
dellammat = diag(dellam,0);

cbarS = inv(sqrt(lammat))*V'*Css1*V*inv(sqrt(lammat));

cbarA = inv(sqrt(Gammamat))*(V'*Css1*delV + delV'*Css1*V + U'*Cas*U)*inv(sqrt(Gammamat)) - 1/2*inv(sqrt(lammat))*V'*Css1*V*lammat^(-3/2)*dellammat - 1/2*lammat^(-3/2)*dellammat*V'*Css1*V*inv(sqrt(lammat));
Cbar = cbarS+cbarA;
%% Calculating Z

[y,omega] = eig(cbarS);

[omega,inds] = sort(diag(omega),1,'descend');

y = y(:, inds);

for k=1:n
    tempind = find(abs(y(:,k))==max(abs(y(:,k))));
    if sign(y(tempind,k)) ~= sign(W_true(tempind,k))
        y(:,k) = y(:,k)*-1;
    end
end

omegamat = diag(omega);

for k=1:n
    delomega(k)=y(:,k)'*cbarA*y(:,k);
    tempvectp = zeros(n,1);
    for p=1:n
        if p==k
            continue
        else
            tempvectp=tempvectp+y(:,p)'*cbarA*y(:,k)/(omega(k)-omega(p))*y(:,p);
        end
    end
    dely(:,k) = tempvectp;
end

delomegamat = diag(delomega);
Phi = omega + delomega';
W = y + dely;

Zs = omegamat^(-1/2)*y'*lammat^(-1/2)*V'*S1;
delZ = omegamat^(-1/2) * y' * lammat^(-1/2) * V' * A1 +...
       omegamat^(-1/2) * y' * lammat^(-1/2) * delV' * S1 +...
       omegamat^(-1/2) * y' * -1/2* lammat^(-3/2)*dellammat * V' * S1 +...
       omegamat^(-1/2) * dely' * lammat^(-1/2) * V' * S1 + ...
       -1/2*omegamat^(-3/2)*delomegamat * y' * lammat^(-1/2) * V' * S1;

Z = Zs + delZ;