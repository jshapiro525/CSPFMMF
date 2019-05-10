clear all;
close all;
clc;
trunc = 1;
imdim = 75;
DATA=[];
name='HR4796a-KL';

for i = [1 20 50 100 999]    
    if i ~= 999
        exten='-speccube.fits';
        filenumb=i;
        numb=num2str(filenumb);
        filename=strcat(name,numb,exten);
    else
        exten='modes-all.fits';
        filename=strcat(name,exten);
    end
    
    data=fitsread(filename,'image');
    if trunc      
        dim=max(size(data));
        c=floor(dim/2);
        w=floor(imdim/2);
        data=data(c-w:c+w,c-w:c+w,:);
    end
    tempdata=sum(data,3);
    DATA=cat(3,DATA,tempdata); 
end 

for i=1:5
    myfig(DATA(:,:,i))
end
