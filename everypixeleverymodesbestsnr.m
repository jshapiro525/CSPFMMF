close all
clear all
clc

load('BPmfmodes')
bploc1=[43 57];
bploc2=[39 18];
planetlocs=[bploc1;bploc2];

% load('HDmfmodes.mat')
% hdloc=[33 58];
% planetlocs=[hdloc];



imdim=75;
nmodes=size(slicemfs,3);
nlams=size(slicemfs,4);
center=[ceil(75/2) ceil(75/2)];

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

modesused = fliplr(1:nmodes);
for pixels = 1:imdim^2
    
    rowloc = rem(pixels,imdim);
    colloc = ceil(pixels/imdim);
    if rowloc<11 || rowloc>imdim-11 || colloc<11 || colloc>imdim-11 || norm([rowloc colloc]-center)<=8
        continue
    end
    
    for k=modesused
        
        if k==nmodes
            [finalsnr(rowloc,colloc),~,~,~]=snrcalc(modes(:,:,k),rowloc,colloc,4,imdim,planetlocs);
            trackingimage=modes(:,:,k);           
            pixelvect(rowloc,colloc,k)=1;
        else
            [tempsnr,tempval,tempdev,tempm]=snrcalc(trackingimage+modes(:,:,k),rowloc,colloc,4,imdim,planetlocs);
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
            if rowloc==33 && colloc == 58 && k == modesused(end)
                thisspot=trackingimage;
            end
          
        end
  
    end
   
end

toc
myfig(finalsnr)
colorbar

myfig(thisspot)

thisspot=thisspot-min(thisspot);
thisspot=thisspot/max(max(thisspot));
myfig(thisspot)


% finalsnr(43,57)
% fdm=[finalval(43,57) finaldev(43,57) finalm(43,57)]

finalsnr(33,58)
fdm=[finalval(33,58) finaldev(33,58) finalm(33,58)]
