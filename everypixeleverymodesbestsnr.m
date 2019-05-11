close all
clear all
clc

load('allmodecspfmmfdata.mat')

imdim=75;
nmodes=size(slicemfs,3);
nlams=size(slicemfs,4);

goodmodes=[1:6 8 9 11 14 15 17 21 23 24 26:29 31 34 35 36];
% goodmodes=[1 2];
goodmodes= 38-goodmodes;

modes=zeros(imdim,imdim,nmodes);

for i = 1:nmodes
    for j = 1:nlams
        modes(:,:,i)=modes(:,:,i)+slicemfs(:,:,i,j);
    end
end


finalsnr=zeros(imdim);

pixelvect=zeros(imdim,imdim,nmodes);


tic

for pixels = 1:imdim^2
    
    rowloc = rem(pixels,imdim);
    colloc = ceil(pixels/imdim);
    if rowloc<10 || rowloc>imdim-10 || colloc<10 || colloc>imdim-10
        continue
    end
    
    for k=fliplr(36:37)
        
        if k==nmodes
            finalsnr(rowloc,colloc)=snrcalc(modes(:,:,k),rowloc,colloc,4,imdim);
            trackingimage=modes(:,:,k);
            pixelvect(rowloc,colloc,k)=1;
        else
            tempsnr=snrcalc(trackingimage+modes(:,:,k),rowloc,colloc,4,imdim);
            if tempsnr>finalsnr(rowloc,colloc)
                finalsnr(rowloc,colloc)=tempsnr;
                pixelvect(rowloc,colloc,k)=1;
                trackingimage=trackingimage+modes(:,:,k);
            else
                pixelvect(rowloc,colloc,k)=0;
            end
          
        end
  
    end
   
end

toc
myfig(finalsnr)
colorbar
