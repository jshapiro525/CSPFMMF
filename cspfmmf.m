close all;
clear all;
clc;

%% Data Creation

epsilon = 0.008; %.008 is used for data analysis
modes = 2; % time scales linearly
% figs = 0;
% 
% load('betpicandinjections75_1.mat')
% temploc = [43 56];

load('hd14706injections75.mat')
temploc = [33 58];

injected = A*epsilon;

% for i=1:length(A(1,1,:));
%     b(:,:,i)=imrotate(injected(:,:,i),90);
% end
totalpert = noise + injected;%+b;

% nlams = 1;
% startlam = 1;
% [totalpert, noise2, injected, lams2] = lamtrunc(totalpert, noise, injected, lams, parangs, nlams, startlam);

a = stackedsnr(totalpert,parangs,lams);

%% Doing CSP 

[CSP,Z_actual,Phi] = splitcspfunc(totalpert, parangs, lams, modes);
% [slicemfs] = truefmmf(Z_actual, noise, A, epsilon, parangs, lams, modes, temploc);

%only one FM way:
figs=0;
[Z, Z_true, Zs, delZ, planetstack, fmstack] = splitFMeveryslice(totalpert, noise, injected, parangs, lams, modes, figs);
for i = 1:length(lams)
     [slicemfs(:,:,:,i)] = fmmf(delZ(:,:,i), temploc, Z_actual(:,:,i), 21, modes);
     fprintf(horzcat(['MF is ' num2str(i/length(lams)*100) ' percent done\n'])) 
end


summfs=zeros(size(slicemfs,1),size(slicemfs,1));
for j = 1:length(lams)
    for i = 1:modes
%         myfig(slicemfs(:,:,i,j))
%         title(horzcat(['mode ' num2str(i) ' wavelength ' num2str(lams(j))]))
        summfs=summfs+slicemfs(:,:,i,j);
    end
end

allmfs = summfs;

for i=1:length(lams)
    allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
end
myfig(summfs)
% fitswrite(allmfs,'allFMMFmap.fits')