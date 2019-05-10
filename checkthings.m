load('HR 4796A.mat')
load('C:/Data/imagedata.mat')
imgs=image_data;
nlams=size(imgs,1);
npars=size(imgs,2);

for lam=1:nlams
    for par = 1:npars
        [a,nshifts]=shiftdim(image_data(lam,par,:,:));
        cspfull(:,:,par,lam)=a;
    end
end

myfig(cspfull(:,:,1,1))
b=sum(cspfull,4);

for k=1:npars
    d(:,:,k) = imrotate(b(:,:,k),parangs(k),'bicubic','crop');
    e(:,:,k) = imrotate(cspfull(:,:,k,1),parangs(k),'bicubic','crop');
end 
c=sum(d,3);
myfig(b(:,:,1))
myfig(c)
myfig(sum(e,3))
eh= b(:,:,1);


for i=1:k
    myfig(b(:,:,i))
end



testsums(:,:,1)=sum(b(:,:,1:4),3);
testsums(:,:,2)=sum(b(:,:,9:16),3);
testsums(:,:,3)=sum(b(:,:,27:29),3);
testsums(:,:,4)=-1*sum(b(:,:,30:37),3);
testsums(:,:,5)=sum(b(:,:,34:37),3);

close all
count=1;
for i=[1 13 28]
    mypars=[1 15 25];
    normed(:,:,count)=b(:,:,i)/norm(b(:,:,i));
    derot(:,:,count) = imrotate(normed(:,:,count),parangs(mypars(count))-parangs(1),'bicubic','crop');
    myfig(normed(:,:,count))
    myfig(derot(:,:,count))
    count=count+1;
end

summed=sum(derot,3);
myfig(summed);
summedall=summed-normed(:,:,2)-normed(:,:,3);
myfig(summedall);

brah(:,:,1)=imrotate(summed,50,'bicubic','crop');
brah(:,:,2)=imrotate(summedall,50,'bicubic','crop');

fitswrite(brah,'C:/Data/finalcspresults.fits')