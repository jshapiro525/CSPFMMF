close all;
clear all;
clc;

%% Data Creation

load('HR 4796A.mat')
DATA=DATA(2:280,2:280,:);
%a = stackedsnr(DATA,parangs,lams);

%% Doing CSP 
total=DATA;
[CSP(:,:,:,:),Z_actual(:,:,:)] = splitcspfunc(total, parangs, lams, 3);


cspfull=CSP;
close all

myfig(cspfull(:,:,1,1))
b=sum(cspfull,4);

for k=1:length(parangs)
    d(:,:,k) = imrotate(b(:,:,k),parangs(k),'bicubic','crop');
    e(:,:,k) = imrotate(cspfull(:,:,k,1),parangs(k),'bicubic','crop');
end 
c=sum(d,3);
myfig(b(:,:,1))
myfig(c)
myfig(sum(e,3))
eh= b(:,:,1);

eh=imrotate(eh,270,'bicubic','crop');
myfig(eh)%



%  for lams=1:length(Z_actual(1,1,:,1))
%      for k=1:length(Z_actual(:,1,1,1))
%         Zfull(:,:,k,lams)=[reshape(Z_actual(k,:,lams,1),imdim,imdim) reshape(Z_actual(k,:,lams,2),imdim,imdim) reshape(Z_actual(k,:,lams,3),imdim,imdim);...
%                            reshape(Z_actual(k,:,lams,4),imdim,imdim) reshape(Z_actual(k,:,lams,5),imdim,imdim) reshape(Z_actual(k,:,lams,6),imdim,imdim);...
%                            reshape(Z_actual(k,:,lams,7),imdim,imdim) reshape(Z_actual(k,:,lams,8),imdim,imdim) reshape(Z_actual(k,:,lams,9),imdim,imdim)];
%      end
%      lams/37
%  end


