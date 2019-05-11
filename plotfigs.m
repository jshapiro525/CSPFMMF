load('allmodecspfmmfdata.mat')

imdim=75;
nmodes=size(slicemfs,3);
nlams=size(slicemfs,4);

goodmodes=[1:6 8 9 11 14 15 17 21 23 24 26:29 31 34 35 36];
% goodmodes=[1 2];
goodmodes= 38-goodmodes;
nmodes=length(goodmodes);

modes=zeros(imdim,imdim,nmodes);

for i = 1:nmodes
    for j = 1:nlams
        modes(:,:,i)=modes(:,:,i)+slicemfs(:,:,i,j);
    end
end

tempentry=[];
for i=1:37
    if ismember(i,goodmodes)
        myfig(modes(:,:,i))
        title(horzcat(['mode ' num2str(i)]))
        colorbar
        
        tempentry=[tempentry max(max(modes(:,:,i)))];
    end
    
end

figure
plot(fliplr(goodmodes),tempentry)