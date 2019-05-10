function [Z,Phi,CSP] = splitcspfunc(Total, parangs, lams, modes)

imdim = size(Total(:,1,1),1); %Side lenth of square image
p=imdim^2; %total number of pixels in an image
n = length(Total(1,1,:)); %total number of images
reflam = lams(ceil(length(lams)/2)); %finding a reference wavelength to scale to

[rep_parangs, rep_lams] = repparlam(parangs,lams); % creates long strings for parangs and lambdas (one for each image)

for thislam=1:length(lams)  % scaling and rotating for each image
    totalperlam = Total(:,:,rep_lams==lams(thislam));
    clear X1
    clear X2
    for thisparang = 1:length(parangs)
        total2(:,:,thisparang)=imrotate(totalperlam(:,:,thisparang),parangs(thisparang),'bicubic','crop');
        X2(:,thisparang)=reshape(total2(:,:,thisparang),p,1);
        X2(:,thisparang)= X2(:,thisparang)-mean(X2(:,thisparang));
        X1(:,thisparang)=reshape(totalperlam(:,:,thisparang),p,1);
        X1(:,thisparang)= X1(:,thisparang)-mean(X1(:,thisparang));
    end
    
    X1=X1';
    X2=X2';

    C1 = X1*X1'*1/(p-1);
    C2 = X2*X2'*1/(p-1);

    Cplus = C1+C2;

    [U,Gamma] = eig(Cplus);
    [Gamma,inds] = sort(diag(Gamma),1,'descend');
    Gamma = diag(Gamma,0);
    U = U(:, inds);

    P = inv(Gamma)^(.5)*U';
    Cbar = P*C1*P';

    [W,Phi] = eig(Cbar);
    [Phi,indsh] = sort(diag(Phi),1,'descend');
    Phi = diag(Phi,0);
    W = W(:, indsh);
    
    Z(:,:,thislam) = sqrt(inv(Phi))*W'*P*X1;
%% ignored
    
    CSP = 0; %remove this for projection results.

    %The next part is for projecting results abck onto original images, but
    %we don't do that.    

%     K = length(parangs)+1-modes;
%     ZK=Z(K:end,:,thislam);
%     temp=(ZK'*ZK);
%     


%     for thisparang=1:length(parangs)
%         sig=temp*X1(thisparang,:)';
%         CSP(:,:,thisparang,thislam)=reshape(sig',imdim,imdim);
%     end
    %fprintf(horzcat(['CSP is ' num2str(thislam/length(lams)*100) ' percent done\n']))
end

end
