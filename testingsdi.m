close all;
clear all;
clc;


%% Data Creation

% Fake Wavelengths and Rotations
 lams = 1:.1:1.9;
% lams = 1;
parangs = [1:2:20];
% parangs = 1;

epsilon = 0.03; % Injected Signal Strength
% imdim = 71; % sets image size to imdim x imdim

% Simple Data
% [total, injected]=generatefakedatalam(lams,parangs,imdim,.05,17.5,epsilon); % "total" represents your final science images, "injected" is just the injected signals
% 

imdim = 115;
load('betpicandinjections.mat')
injected = A*epsilon;
total = noise + injected;


% Beta Pic data
%[total, parangs, lams, centers, injected] = generateinjectedsamer(imdim,17,1,sigstr,'BetPic');

% Beta Pic just SDI:
% total = total(:,:,1:length(lams));
% parangs = 1;

% Beta Pic just ADI:
% for i = 1:length(parangs)
%     aditotal(:,:,i) = total(:,:,1 + length(lams)*(i-1));
% end
% total = aditotal;
% lams = 1;

disp('Data Created')

a = stackedsnr(total,parangs,lams);

%% Doing CSP 

[CSP, eigenvalues, Z, K, x1, x2] = cspfunctionallrot(total, parangs, lams, imdim, 50);

%% Different figures for analysis

figure()
plot(eigenvalues)

myfig(injected(:,:,1))

% for i=1:length(x1(1,1,:))
%     myfig(x1(:,:,i))
%     title(horzcat(['SDI Cropped Wavelength: ' num2str(lams(i))]))
% end

% for i=1:length(CSP(1,1,:))
%     myfig(CSP(:,:,i))
%     title(horzcat(['SDI Image ' num2str(i)]))
% end

% for i=1:length(CSP(1,1,:))
%     myfig(reshape(Z(i,:),imdim,imdim))
%     title(horzcat(['mode ' num2str(i)]))
% end
  
finalstack = realign(CSP,parangs,lams);

myfig(finalstack)
title('Result')