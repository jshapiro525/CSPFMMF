function [Z_true, Phi_true, Z, Zs, delZ] = splitFM(noise, injected, modes, parangs, figs)

    totalpert = noise + injected;
    imdim = length(totalpert(:,1,1));
    pix = imdim^2; %total number of pixels in an image
    n = length(totalpert(1,1,:)); %total number of images

    %% Getting actual cov matrix and Z from perturbed data

    X1 = [];
    X2 = [];

    for i=1:n  % scaling and rotating for each image

        total2(:,:,i)=imrotate(totalpert(:,:,i),parangs(i),'bicubic','crop');
        X2(:,i)=reshape(total2(:,:,i),pix,1);
        X1(:,i)=reshape(totalpert(:,:,i),pix,1);
        X2(:,i)= X2(:,i)-mean(X2(:,i));
        X1(:,i)= X1(:,i)-mean(X1(:,i));
    end

    X1=X1';
    X2=X2';

    C1 = X1*X1'/(pix-1);
    C2 = X2*X2'/(pix-1);%*1/trace(X2*X2');

    Cplus = C1+C2;

    [U_true,Gamma_true] = eig(Cplus);
    [Gamma_true,inds] = sort(diag(Gamma_true),1,'descend');
    Gamma_true = diag(Gamma_true,0);
    U_true = U_true(:, inds);

%     for k=1:n
%         if U_true(1,k) < 0
%             U_true(:,k) = U_true(:,k)*-1;
%         end
%     end

    P_true = inv(Gamma_true)^(0.5)*U_true';
    Cbar_true = P_true*C1*P_true';

    [W_true,Phi_true] = eig(Cbar_true);
    [Phi_true,indsh] = sort(diag(Phi_true),1,'descend');
    W_true = W_true(:, indsh);
% 
%     for k=1:n
%         if W_true(1,k) < 0
%             W_true(:,k) = W_true(:,k)*-1;
%         end
%     end

    Z_true = sqrt(inv(diag(Phi_true)))*W_true'*P_true*X1;

    %% unperturbed Cov, eigenvals, and Z

    for i=1:n  % scaling and rotating for each image

        total2(:,:,i)=imrotate(noise(:,:,i),parangs(i),'bicubic','crop');
        S2(:,i)=reshape(total2(:,:,i),pix,1);
        S1(:,i)=reshape(noise(:,:,i),pix,1);
        S2(:,i)= S2(:,i)-mean(S2(:,i));
        S1(:,i)= S1(:,i)-mean(S1(:,i));

    end

    S1=S1';
    S2=S2';

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

    %% Signals and their perturbations

    for i=1:n  % scaling and rotating for each image

        total2(:,:,i)=imrotate(injected(:,:,i),parangs(i),'bicubic','crop');
        A2(:,i)=reshape(total2(:,:,i),pix,1);
        A1(:,i)=reshape(injected(:,:,i),pix,1);
        A2(:,i)= A2(:,i)-mean(A2(:,i));
        A1(:,i)= A1(:,i)-mean(A1(:,i));
    end

    A1=A1';
    A2=A2';

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

    
    %% figures
    if figs
        Gamma_true = diag(Gamma_true);
        
%         figure()
%         plot(1:n,Gamma_true,'r',1:n,Gamma,'b')
%         legend('True Gamma','Estimated Gamma')
%         
%         figure()
%         plot(abs(Gamma_true-Gamma)./Gamma_true*100)
%         title('Percent Gamma (Cplus Eigenvalue) Difference')

        pdifU= abs((U_true-U)./U_true)*100;
        log_graph_element_distro(pdifU,'Percent U (Cplus Eigenvector) Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

        pdifV= abs((U_true-V)./U_true)*100;
        myfig(pdifV)
        log_graph_element_distro(pdifV,'Percent V (Cplus Eigenvector) Difference','Percent Error','Number of entries in U','Distribution of Percent Error')
        
%         perdifcbar = abs((Cbar-Cbar_true)./Cbar_true)*100; 
%         log_graph_element_distro(perdifcbar,'Percent FM Cbar Difference','Percent Error','Number of entries in U','Distribution of Percent Error')
%         
%         perdifcbarun = abs((Cbar_un-Cbar_true)./Cbar_true)*100; 
%         log_graph_element_distro(perdifcbarun,'Percent unpert FM Cbar Difference','Percent Error','Number of entries in U','Distribution of Percent Error')

%         figure()
%         plot(1:n,Phi_true,'r',1:n,Phi,'b')
%         legend('True Phi','Estimated Phi')
%         
%         figure()
%         plot(abs((Phi_true-Phi)./Phi_true*100))
%         title('Percent Phi (Cbar Eigenvalue) difference')
%         
%         figure()
%         plot(abs((Phi_true-omega)./Phi_true*100))
%         title('Percent omega (Cbar Eigenvalue) difference')
        
        
         pdifW = abs((W_true-W)./W_true)*100;
%         myfig(pdifW)
%         title('Percent w (Cbar Eigenvector) Difference')  
%         
         pdifY = abs((W_true-y)./W_true)*100;
%         myfig(pdifY)
%         title('Percent Y (Cbar Eigenvector) Difference')
        
        log_graph_element_distro(pdifW,'Distribution of W (Cbar Eigenvector) percent error','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');
%         log_graph_element_distro(pdifY,'Distribution of Y (Cbar Eigenvector) percent error','Percent Error','Number of entries in the eigenvectors','Distribution of Percent Error');

        pdifZ = abs((Z_true-Z)./Z)*100;
        pdifZun = abs((Z_true-Zs)./Z_true)*100;
        log_graph_element_distro(pdifZun,'Difference in perturbed projection matrices to unperturbed model','Percent Error','Number of entries in Z','Distribution of Percent Error');

        log_graph_element_distro(pdifZ,'Difference in perturbed projection matrices','Percent Error','Number of entries in Z','Distribution of Percent Error');
    end
end
