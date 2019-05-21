function [snr,val,dev,m] = snrcalc(Image,i,j,annulussize,imdim,planetlocs)

loc=[i j];
val=Image(loc(1),loc(2));
center=[ceil(imdim/2) ceil(imdim/2)];
    
for i=1:imdim
    for j=1:imdim
        if norm([i j]-center)<=10 || norm([i j]-loc)<=5
            Image(i,j)=nan;
        end
    end
end

cropped=Image(11:imdim-11,11:imdim-11);

[rows, cols]=size(cropped);
[rr,cc] = meshgrid(1:size(cropped));

baseimage=cropped;
r= norm(center-loc);
mask = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)>= (r-annulussize/2);
mask2 = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)<= (r+annulussize/2);
mask_annulus = mask&mask2;

maskall=boolean(ones(size(baseimage,1)));
for i =1: size(planetlocs,1)
    for k=1:size(baseimage,1)
        for j=1:size(baseimage,1)
            if norm(([k j])+11-planetlocs(i,:))<=5
                maskall(k,j)=0;
            end
        end
    end
    maskall = maskall&mask_annulus;
end
    
    
masked=(baseimage+eps).*maskall;
masked(masked==0)=nan;
vected=reshape(masked-eps,length(cropped(:,1))^2,1);
vected=vected(~isnan(vected));
m=mean(vected);
dev=std(vected);
snr = (val-m)./dev;

end