close all;
clear all;
clc;

%% Data Creation

load('HR 4796A.mat')
DATA=DATA(2:280,2:280,:);
%a = stackedsnr(DATA,parangs,lams);

%% Doing CSP 

for i=1:9
    if i==1
        dim1=[1:93];
        dim2=dim1;
        total=DATA(dim1,dim2,:);
    elseif i==2
        dim1=[1:93];
        dim2=[94:186];
        total=DATA(dim1,dim2,:);
    elseif i==3
        dim1=[1:93];
        dim2=[187:279];
        total=DATA(dim1,dim2,:);
    elseif i==4
        dim1=[94:186];
        dim2=[1:93];
        total=DATA(dim1,dim2,:);
    elseif i==5
        dim1=[94:186];
        dim2=[94:186];
        total=DATA(dim1,dim2,:);
    elseif i==6
        dim1=[94:186];
        dim2=[187:279];
        total=DATA(dim1,dim2,:);        
    elseif i==7
        dim1=[187:279];
        dim2=[1:93];
        total=DATA(dim1,dim2,:);
    elseif i==8
        dim1=[187:279];
        dim2=[94:186];
        total=DATA(dim1,dim2,:);
    else
        dim1=[187:279];
        dim2=[187:279];
        total=DATA(dim1,dim2,:);            
    end
        
    [CSP(:,:,:,:,i),Z_actual(:,:,:,i),Phi] = splitcspfunc(total, parangs, lams, 2);
end

total=DATA;
% [CSP(:,:,:,:,i),Z_actual(:,:,:,i)] = splitcspfuncparl(total, parangs, lams, 3);

%to start here, execute following command:
load('hr4796a after.mat')

cspfull = [CSP(:,:,:,:,1) CSP(:,:,:,:,2) CSP(:,:,:,:,3);...
            CSP(:,:,:,:,4) CSP(:,:,:,:,5) CSP(:,:,:,:,6);...
            CSP(:,:,:,:,7) CSP(:,:,:,:,8) CSP(:,:,:,:,9)];
 
% cspfull=CSP;
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
%  
%  for lams=1:length(Z_actual(1,1,:,1))
%      for k=1:length(Z_actual(:,1,1,1))
%         Zfull(:,:,k,lams)=[reshape(Z_actual(k,:,lams,1),imdim,imdim) reshape(Z_actual(k,:,lams,2),imdim,imdim) reshape(Z_actual(k,:,lams,3),imdim,imdim);...
%                            reshape(Z_actual(k,:,lams,4),imdim,imdim) reshape(Z_actual(k,:,lams,5),imdim,imdim) reshape(Z_actual(k,:,lams,6),imdim,imdim);...
%                            reshape(Z_actual(k,:,lams,7),imdim,imdim) reshape(Z_actual(k,:,lams,8),imdim,imdim) reshape(Z_actual(k,:,lams,9),imdim,imdim)];
%      end
%      lams/37
%  end

eh= b(:,:,1);

eh=imrotate(eh,270,'bicubic','crop');
myfig(eh)
