close all;
clear all;
clc;

tic

center=[38 38];

for k=1:2
    if k==1
        load('Allresultsworkspace.mat')
        load('betapicworkspaceres.mat')

    end
    
    if k==2
        [summfs1, summfs]=klipfmmfresultssaved(0);
    end
        
    
    
%     summfs = summfs-min(min(summfs));
%     summfs = summfs/max(max(summfs))*100;
% 
%     summfs1 = summfs1-min(min(summfs1));
%     summfs1 = summfs1/max(max(summfs1))*100;
    
   
%      myfig(summfs1)
%      myfig(summfs)
   
    if k==1
        bploc1=[43 57];
        bploc2=[39 18];
        hdloc=[33 58];
        vals(1:3)=[summfs1(bploc1(1),bploc1(2)) summfs1(bploc2(1),bploc2(2)) summfs(hdloc(1),hdloc(2))];
    else
        bploc1=[46 55];
        bploc2=[40 21];
        hdloc=[34 57];
        vals(4:6)=[summfs1(bploc1(1),bploc1(2)) summfs1(bploc2(1),bploc2(2)) summfs(hdloc(1),hdloc(2))];
        bploc2=[39 19];
    end
    sizen=75;
    for i= 11:sizen-11
        for j=11:sizen-11
            if k==1
                if norm([i j]-center)<=5
                    summfs1(i,j)=nan;
                    summfs(i,j)=nan;
                end
            else
                if norm([i j]-center)<=8
                    summfs1(i,j)=nan;
                end
                if norm([i j]-fliplr(center))<=8
                    summfs(i,j)=nan;
                end
            end
            if norm([i j]-bploc1)<=5 || norm([i j]-bploc2)<=5
                summfs1(i,j)=nan;
            end
            if norm([i j]-hdloc)<=5
                summfs(i,j)=nan;
            end
        end
    end
    
    betpic=summfs1(11:sizen-11,11:sizen-11);
    hd=summfs(11:sizen-11,11:sizen-11);
% 
%      myfig(betpic)
%      myfig(hd)


    annulussize = 4;

    [rows, cols]=size(betpic);
    [rr,cc] = meshgrid(1:size(betpic));

    for i=1:3
        baseimage=betpic;
        if i==1
            r= norm(center-bploc1);
        elseif i==2
            r= norm(center-bploc2);
        else
            r= norm(center-hdloc);
            if k==2
                r= norm(fliplr(center)-hdloc);
            end
            baseimage=hd;
        end
        mask = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)>= (r-annulussize/2);
        mask2 = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)<= (r+annulussize/2);
        mask_annulus = mask&mask2;
        masked=(baseimage+.0000000001).*mask_annulus;
        masked(masked==0)=nan;
          myfig(masked);
        vected=reshape(masked,length(betpic(:,1))^2,1);
        vected=vected(~isnan(vected));
       
        content=0;
        means((k-1)*3+i)=mean(vected);
        mvected=vected-mean(vected);
        mvals((k-1)*3+i)=vals((k-1)*3+i)-mean(vected);
        mstd=std(mvected);
        
        for m=1:length(vected)
            if abs(mvected(m))<= mstd
                content=content+1;
            end
        end
        
        standarddev((k-1)*3+i)=mstd;
        skew((k-1)*3+i)=skewness(vected);
        kurt((k-1)*3+i)=kurtosis(vected);
        percent((k-1)*3+i)=content/length(vected);
        
        %[cdf((k-1)*3+i),esnr1((k-1)*3+i),lambda((k-1)*3+i)] = probmapper(vals((k-1)*3+i),vected);
        
        disp(horzcat([num2str((k-1)*3+i) ' of 6 Finshed!' newline]))
           
     end
%     
%    
% %     ADres
% %     p
end
% type
snr = (vals-means)./standarddev
cdf
esnr1


toc



 