close all;
clear all;
clc;

%% Data Creation

epsilon = 0.05;
modes = 2; 

%% HD 
tic

load('FinalHDWorkSpace.mat')
temploc = [46 71];

injected = A*epsilon;
noise = totalpert - injected;

a = stackedsnr(totalpert,parangs,lams);

[injtrunc, noisetrunc, ttrunc] = splituppic(injected, noise, totalpert);
for i=1:9
    [CSP,Z_actual(:,:,:,i),Phi] = splitcspfunc(ttrunc(:,:,:,i), parangs, lams, modes);
    fprintf(horzcat(['CSP is ' num2str(i/9*100) ' percent done\n']))
end

[slicemfs] = segmentedtruefmmf(Z_actual, noise, injected, totalpert, parangs, lams, modes, temploc);

summfs=zeros(size(slicemfs,1),size(slicemfs,1));
for j = 1:length(lams)
    for i = 1:modes
        summfs=summfs+slicemfs(:,:,i,j);        
    end
end

allmfs = summfs;

for i=1:length(lams)
    allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
end
myfig(summfs)
toc
% fitswrite(allmfs,'allFMMFmap.fits')

%% BP

% load('FinalBPWorkSpace.mat')
% temploc=[56 69];