close all;
clear all;
clc;

%% Data Creation

epsilon = 0.03;
modes = 2;
figs = 1;

% lams = 1.2:(1.8-1.2)/36:1.8;
% lams = 1;
% parangs = [1:31/36:32];
% imdim = 71;
% [totalpert, injected, noise]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon);

imdim = 75;
%load('hd14706injections75.mat')
load('betpicandinjections75_1.mat')
injected = A*epsilon;
totalpert = noise + injected;
npar = length(parangs);

nlams = 2;
lams=lams(1:nlams);
for i=1:npar
    tert(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) =totalpert(:,:,(i-1)*37+[1:nlams]); 
    nois(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = noise(:,:,(i-1)*37+[1:nlams]);
    inj(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = injected(:,:,(i-1)*37+[1:nlams]);
end
totalpert=tert;
noise=nois;
injected=inj;
% 
% nlams=1;
% lams=lams(1:nlams);
% for i=1:npar
%     tert(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) =totalpert(:,:,(i-1)*37+[15:15+nlams-1]); 
%     nois(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = noise(:,:,(i-1)*37+[15:15+nlams-1]);
%     inj(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = injected(:,:,(i-1)*37+[15:15+nlams-1]);
% end
% totalpert=tert;
% noise=nois;
% injected=inj;


% a = stackedsnr(totalpert,parangs,lams);

%% Doing CSP 

[CSP,Z_actual,Phi] = splitcspfunc(totalpert, parangs, lams, modes);
[Z, Z_true, Zs, delZ, planetstack, fmstack] = splitFMeveryslice(totalpert, noise, injected,parangs, lams, modes, figs);

% log_graph_element_distro(abs((Z_actual-Z_true)./Z_actual)*100,'thistitle','number of','percent diff','does Z match?')

% %% Different figures for analysis
% 
% myfig(injected(:,:,1))
% title('injected location')
% % 
for i=1:length(lams)
    cspstack(:,:,i) = realignFM(CSP(:,:,:,i),parangs,lams);
    %myfig(finalstack(:,:,i))
    %title(horzcat(['stack for wavelength ' num2str(lams(i))]))

%     myfig(cspstack(:,:,i))
%     title(horzcat(['stack for wavelength ' num2str(lams(i))]))
% % 
%     myfig(fmstack(:,:,i))
%     title(horzcat(['FM stack for wavelength ' num2str(lams(i))]))
% 
%     myfig(planetstack(:,:,i))
%     title(horzcat(['Planet FM stack for wavelength ' num2str(lams(i))]))
end
% % % 
myfig(sum(cspstack,3))
%title('All Summed')
myfig(sum(fmstack,3))
%title('FM Sum')

betapiccomparisonplotter;

myfig(sum(planetstack,3))
%title('planetFM Sum')