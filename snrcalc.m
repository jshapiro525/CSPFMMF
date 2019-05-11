function [snr] = snrcalc(Image,i,j,annulussize,imdim)

loc=[i j];
val=Image(loc(1),loc(2));
center=[ceil(imdim/2) ceil(imdim/2)];
    
for i= center(1)-7:center(1)+7
    for j=center(2)-7:center(2)+7
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
masked=(baseimage+eps).*mask_annulus;
masked(masked==0)=nan;
vected=reshape(masked-eps,length(cropped(:,1))^2,1);
vected=vected(~isnan(vected));
snr = (val-mean(vected))./std(vected);

end