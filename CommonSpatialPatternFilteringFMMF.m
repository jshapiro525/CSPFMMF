close all;
clear all;
clc;

%% Beta Pic
tic
% Data Creation
epsilon = 0.008; %.008 is used for data analysis
modes = 2; 

load('betpicandinjections75_1.mat')
temploc = [43 56];

% load('hd14706injections75.mat')
% temploc = [33 58];

injected = A*epsilon;
totalpert = noise + injected;

% Doing CSP and FM

% [Z_actual,Phi,~] = splitcspfunc(totalpert, parangs, lams, modes); %Done
% in splitfmeveryslice

figs=0;
[Z_true, Phi_true, Z, Zs, delZ, planetstack, fmstack] = splitFMeveryslice(totalpert, noise, injected, parangs, lams, modes, figs);
for i = 1:length(lams)
     [slicemfs(:,:,:,i)] = fmmf(delZ(:,:,i), temploc, Z_true(:,:,i), 21, modes);
     fprintf(horzcat(['MF is ' num2str(i/length(lams)*100) ' percent done\n'])) 
end


% Combining data
summfs1=zeros(size(slicemfs,1),size(slicemfs,1));
for j = 1:length(lams)
    for i = 1:modes
%         myfig(slicemfs(:,:,i,j))
%         title(horzcat(['mode ' num2str(i) ' wavelength ' num2str(lams(j))]))
        summfs1=summfs1+slicemfs(:,:,i,j);
    end
end

allmfs = summfs1;

for i=1:length(lams)
    allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
end

% Plotting Figures

myfig(injected(temploc(1)-5:temploc(1)+5,temploc(2)-5:temploc(2)+5,1))
myfig(injected(:,:,1))
a = stackedsnr(totalpert,parangs,lams);

myfigcbar(reshape(Z_true(end-1,:,12)+Z_true(end,:,12),75,75))
myfigcbar(summfs1)
save('betapicworkspaceres','summfs1')

%% HD 14706
clear all;

% Data Creation

load('betapicworkspaceres')

epsilon = 0.008; %.008 is used for data analysis
modes = 2; 

load('hd14706injections75.mat')
temploc = [33 58];

injected = A*epsilon;
totalpert = noise + injected;

% Doing CSP and FM

[CSP,Z_actual,Phi] = splitcspfunc(totalpert, parangs, lams, modes);

figs=0;
[Z_true, Phi_true, Z, Zs, delZ, planetstack, fmstack] = splitFMeveryslice(totalpert, noise, injected, parangs, lams, modes, figs);
for i = 1:length(lams)
     [slicemfs(:,:,:,i)] = fmmf(delZ(:,:,i), temploc, Z_true(:,:,i), 21, modes);
     fprintf(horzcat(['MF is ' num2str(i/length(lams)*100) ' percent done\n'])) 
end


% Combining data
summfs2=zeros(size(slicemfs,1),size(slicemfs,1));
for j = 1:length(lams)
    for i = 1:modes
%         myfig(slicemfs(:,:,i,j))
%         title(horzcat(['mode ' num2str(i) ' wavelength ' num2str(lams(j))]))
        summfs2=summfs2+slicemfs(:,:,i,j);
    end
end

allmfs = summfs2;

for i=1:length(lams)
    allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
end

% Plotting Figures

a = stackedsnr(totalpert,parangs,lams);

myfigcbar(summfs2)

pdifZ = abs((Z_true(:,:,12)-Z(:,:,12))./Z(:,:,12))*100;
log_graph_element_distro(pdifZ,'','Percent Error','Number of entries in Z','Distribution of Percent Error');

save('hdworkspaceres','summfs2')
klipfmmfresults()
toc
