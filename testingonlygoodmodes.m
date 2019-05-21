close all;
clear all;
clc;

%% Beta Pic
tic
% Data Creation
epsilon = 0.008; %.008 is used for data analysis

load('betpicandinjections75_1.mat')
temploc = [43 57];

% load('hd14706injections75.mat')
% temploc = [33 58];

injected = A*epsilon;
totalpert = noise + injected;

goodmodes=[1:6 8 9 11 14 15 17 21 23 24 26:29 31 34 35 36];
% goodmodes=[1 2];
goodmodes= 38-goodmodes;
nmodes=length(goodmodes);

figs=0;
[Z_true, Phi_true, Z, Zs, delZ] = splitFMeveryslice(totalpert, noise, injected, parangs, lams, goodmodes, figs);


for i = 1:length(lams)
     [slicemfs(:,:,:,i)] = fmmf(delZ(:,:,i), temploc, Z_true(:,:,i), 21, goodmodes);
     fprintf(horzcat(['MF is ' num2str(i/length(lams)*100) ' percent done\n'])) 
end


% Combining data
summfsall=zeros(size(slicemfs,1),size(slicemfs,1));
for j = 1:length(lams)
    for i = 1:nmodes
        summfsall=summfsall+slicemfs(:,:,i,j);
    end
end

allmfs = summfsall;


for i=1:length(lams)
    allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
end
toc

center=[38 38];


finalimage=summfsall;

[i1,i2]=find(finalimage==max(max(finalimage(:,ceil(end/2):end))));
bploc1=[i1 i2];
[i1,i2]=find(finalimage==max(max(finalimage(:,1:ceil(end/2)))));
bploc2=[i1 i2];
%         hdloc=[33 58];
vals=[finalimage(bploc1(1),bploc1(2)) finalimage(bploc2(1),bploc2(2))];% finalimage(hdloc(1),hdloc(2))];

sizen=75;
for i= 11:sizen-11
    for j=11:sizen-11
        if norm([i j]-center)<=5 || norm([i j]-bploc1)<=5 || norm([i j]-bploc2)<=5
            finalimage(i,j)=nan;
        end
%             if norm([i j]-hdloc)<=5
%                 summfs(i,j)=nan;
%             end
    end
end

cropped=finalimage(11:sizen-11,11:sizen-11);


annulussize = 4;

[rows, cols]=size(cropped);
[rr,cc] = meshgrid(1:size(cropped));

for i=1:2
    baseimage=cropped;
    if i==1
        r= norm(center-bploc1);
    elseif i==2
        r= norm(center-bploc2);
    end
    mask = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)>= (r-annulussize/2);
    mask2 = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)<= (r+annulussize/2);
    mask_annulus = mask&mask2;
    masked=(baseimage+.0000000001).*mask_annulus;
    masked(masked==0)=nan;
      %myfig(masked);
    vected=reshape(masked,length(cropped(:,1))^2,1);
    vected=vected(~isnan(vected));

    means(i)=mean(vected);
    standarddev(i)=std(vected);         

 end

% type
snr = (vals-means)./standarddev
% figure
% plot([1:maxnmodes;1:maxnmodes]',snr);
% xlabel('Number of Modes Used')
% ylabel('SNR')
% legend('injected planet','beta pic b')

% save('HDmultimodedata')
snr = (vals)./standarddev
