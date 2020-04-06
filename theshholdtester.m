close all
clear all
clc

%% Loading, plotting original
load('C:/Users/jshap/Documents/pyCSP/testmodeldata.mat')
n=size(models,1);%total number of images
imdim=281;
pix = size(A1,2);


dZfig=nan(n,imdim^2);
dZfig(:,inds) = dZ;
for i=1:n
    dZfigresh(:,:,i)=reshape(dZfig(i,:),imdim,imdim);
end
myfig(dZfigresh(:,:,n))
title('unthresholded delta_z')





%% testing threshold


for thresh = [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5]   
    dZfig=nan(n,imdim^2);
    for i=1:n
        dZstd=std(dZ(i,:));
        dZtemp=dZ(i,:);
        dZtemp(find(abs(dZ(i,:)-median(dZ(i,:)))<thresh*dZstd))=nan;    
        dZfig(i,inds) = dZtemp;
        dZfigthresh(:,:,i)=reshape(dZfig(i,:),imdim,imdim);
    end
%     myfig(dZfigthresh(:,:,n))
%     title(horzcat(['Threshold = ' num2str(thresh) ' delta_z']))
    fitswrite(dZfigthresh,horzcat(['dZ_thresh_' num2str(thresh) '.fits']))
end

