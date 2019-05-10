function [CSP,Z] = splitcspfuncparl(Total, parangs, lams, modes)

imdim = size(Total(:,1,1),1); %Side lenth of square image
p=imdim^2; %total number of pixels in an image
n = length(Total(1,1,:)); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to

[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strings for parangs and lambdas (one for each image)
nlams=length(parangs);
K = length(parangs)+1-modes;

parfor thislam=1:length(lams)  % scaling and rotating for each image
    totalperlam = Total(:,:,rep_lams==lams(thislam));

    X1=[];
    X2=[];
    total2=[];
    for thisparang = 1:nlams
        total2(:,:,thisparang,thislam)=imrotate(totalperlam(:,:,thisparang),parangs(thisparang),'bicubic','crop');
        X2t=reshape(total2(:,:,thisparang,thislam),1,p);
        X2t=X2t-mean(X2t);
        X1t=reshape(totalperlam(:,:,thisparang),1,p);
        X1t= X1t-mean(X1t);
        X2(thisparang,:) = X2t;
        X1(thisparang,:) = X1t;   
    end
    
    C1 = X1*X1'*1/(p-1);
    C2 = X2*X2'*1/(p-1);

    Cplus = C1+C2;

    [U,Gamma] = eig(Cplus);
    [Gamma,inds] = sort(diag(Gamma),1,'descend');
    Gamma = diag(Gamma,0);
    U = U(:, inds);

    P = inv(Gamma)^(.5)*U';
    Cbar = P*C1*P';

    [W,Phi(:,:,thislam)] = eig(Cbar);
    b=Phi(:,:,thislam);
    [b,indsh] = sort(diag(b),1,'descend');
    b = diag(b,0);
    W = W(:, indsh);
    
    Z(:,:,thislam) = sqrt(inv(b))*W'*P*X1;
    d = Z(:,:,thislam);
    ZK = d(K:K+modes-1,:);
    temp=(ZK'*ZK);
    
    g=[];
    for thisparang=1:nlams
        sig=temp*X1(thisparang,:)';
        g(:,:,thisparang,thislam)=reshape(sig',imdim,imdim);
    end
    
    CSP(:,:,:,thislam) = g(:,:,:,thislam);
    
    fprintf(horzcat(['CSP is ' num2str(thislam/length(lams)*100) ' percent done\n']))
end

end
