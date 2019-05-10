close all;
clear all;
clc;

%% Data Creation

epsilon = 0.08;
modes = 2; 

%% HD 

load('FinalBPWorkSpace.mat')
temploc = [56 69]-18;

A=A(19:end-18,19:end-18,:);
totalpert=DATA(281/2+.5-32:281/2+.5+32,281/2+.5-32:281/2+.5+32,:);

%temploc = [46 71];

injected = A*epsilon;
noise = totalpert - injected;

a = stackedsnr(totalpert,parangs,lams);


[CSP,Z_actual,Phi] = splitcspfunc(totalpert, parangs, lams, modes);
[slicemfs] = truefmmf(Z_actual, noise, injected, parangs, lams, modes, temploc);

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