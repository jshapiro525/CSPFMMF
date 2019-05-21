close all
clear all
clc

load('allmodecspfmmfdata.mat')

imdim=75;
nmodes=size(slicemfs,3);
nlams=size(slicemfs,4);

bploc1=[43 57];
bploc2=[39 18];

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

modesused = fliplr(1:37);
for pixels = 1:imdim^2
    
    rowloc = rem(pixels,imdim);
    colloc = ceil(pixels/imdim);
    if rowloc<11 || rowloc>imdim-11 || colloc<11 || colloc>imdim-11
        continue
    end
    
    for k=modesused
        
        if k==nmodes
            [finalsnr(rowloc,colloc),~,~,~]=snrcalc(modes(:,:,k),rowloc,colloc,4,imdim,[bploc1;bploc2]);
            trackingimage=modes(:,:,k);           
            pixelvect(rowloc,colloc,k)=1;
        else
            [tempsnr,tempval,tempdev,tempm]=snrcalc(trackingimage+modes(:,:,k),rowloc,colloc,4,imdim,[bploc1;bploc2]);
            if tempsnr>finalsnr(rowloc,colloc)
                finalsnr(rowloc,colloc)=tempsnr;
                finalval(rowloc,colloc)=tempval;
                finaldev(rowloc,colloc)=tempdev;
                finalm(rowloc,colloc)=tempm;
                pixelvect(rowloc,colloc,k)=1;
                trackingimage=trackingimage+modes(:,:,k);
            else
                pixelvect(rowloc,colloc,k)=0;
            end
            if rowloc==43 && colloc == 57 && k == modesused(end)
                thisspot=trackingimage;
            end
          
        end
  
    end
   
end

toc
myfig(finalsnr)
colorbar

myfig(thisspot)
colorbar

finalsnr(43,57)
fdm=[finalval(43,57) finaldev(43,57) finalm(43,57)]
