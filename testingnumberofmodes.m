close all;
clear all;
clc;

%% Beta Pic
tic
% Data Creation
epsilon = 0.008; %.008 is used for data analysis

% load('betpicandinjections75_1.mat')
% temploc = [43 56];

load('hd14706injections75.mat')
temploc = [33 58];

injected = A*epsilon;
totalpert = noise + injected;

maxnmodes=10;

figs=0;
[Z_true, Phi_true, Z, Zs, delZ, planetstack, fmstack] = splitFMeveryslice(totalpert, noise, injected, parangs, lams, maxnmodes, figs);

for nmodes=1:maxnmodes
    clear slicemfs    
    for i = 1:length(lams)
         [slicemfs(:,:,:,i)] = fmmf(delZ(:,:,i), temploc, Z_true(:,:,i), 21, nmodes);
         fprintf(horzcat(['MF is ' num2str(i/length(lams)*100) ' percent done\n'])) 
    end


    % Combining data
    summfs1=zeros(size(slicemfs,1),size(slicemfs,1));
    for j = 1:length(lams)
        for i = 1:nmodes
            summfs1=summfs1+slicemfs(:,:,i,j);
        end
    end

    allmfs = summfs1;

    for i=1:length(lams)
        allmfs = cat(3,allmfs, slicemfs(:,:,:,i));
    end
    
    summfsall(:,:,nmodes)=summfs1;

    disp(horzcat(['WE ARE THROUGH ' num2str(nmodes) '/' num2str(maxnmodes) ' NMODE SELECTIONS' newline newline newline]))
end
toc

center=[38 38];

for k=1:maxnmodes

    finalimage=summfsall(:,:,k);
         
    [i1,i2]=find(finalimage==max(max(finalimage(:,ceil(end/2):end))));
    bploc1=[i1 i2];
    [i1,i2]=find(finalimage==max(max(finalimage(:,1:ceil(end/2)))));
    bploc2=[i1 i2];
%         hdloc=[33 58];
    vals(k,:)=[finalimage(bploc1(1),bploc1(2)) finalimage(bploc2(1),bploc2(2))];% finalimage(hdloc(1),hdloc(2))];

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
       
        means(k,i)=mean(vected);
        standarddev(k,i)=std(vected);         
        
     end

end
% type
snr = (vals-means)./standarddev;
figure
plot([1:maxnmodes;1:maxnmodes]',snr);
xlabel('Number of Modes Used')
ylabel('SNR')
legend('injected planet','beta pic b')

% save('HDmultimodedata')
snr = (vals)./standarddev
