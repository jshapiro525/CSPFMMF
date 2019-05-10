close all;
clear all;
clc;

import matlab.io.*

imdim=75;
epsilon=.08;
DATA =[];
load('hd14706injections75.mat') % hd
%load('betpicandinjections75_1.mat') % bp
injected = A*epsilon;
totalpert = noise + injected;

        name='S20170806S0'; % HD
%        name='S20151106S0'; % bp
        exten='_spdc_distorcorr.fits';
 for i = 266:303 % hd
     if i~=268 & i~=273 & i~=274 % HD
%for i = 291:328 % bp
%        if i~=304 % bp
        filenumb=i;
        if filenumb<10
            numb=strcat('00',num2str(filenumb));
        elseif filenumb<100
            numb=strcat('0',num2str(filenumb));
        else
            numb=num2str(filenumb);
        end
        filename=strcat(name,numb,exten);
        data=fitsread(filename,'image');
        DATA=cat(3,DATA,data); 

    end
end    

for i=1:size(data,1)
    for j=1:size(data,2)
        for k=1:size(DATA,3)
            if isnan(DATA(i,j,k))
                DATA(i,j,k)=0;
            end
        end
    end
end

DATA = DATA./max(max(max(DATA)));

count=1;
dim=max(size(data,1));
c=floor(dim/2)+4; %hd
%c=floor(dim/2)+6; %bp
w=floor(imdim/2);

for i=1:length(DATA(1,1,:))
    DATA(c-w:c+w,c-w:c+w,i)=totalpert(:,:,i);
    newdata(i,:,:)=DATA(:,:,i);
    if rem(i,100) == 0 | i==length(DATA(1,1,:))
        i/length(DATA(1,1,:))
    end
end