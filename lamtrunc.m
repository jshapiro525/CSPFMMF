function [totalpert, noise, injected, lams] = lamtrunc(totalpert, noise, injected, lams, parangs, nlams, startlam)
npar=length(parangs);
lams=lams(startlam:nlams+startlam-1);
for i=1:npar
    tert(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) =totalpert(:,:,(i-1)*37+[startlam:nlams+startlam-1]); 
    nois(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = noise(:,:,(i-1)*37+[startlam:nlams+startlam-1]);
    inj(:,:,nlams*(i-1)+1:nlams*(i-1)+nlams) = injected(:,:,(i-1)*37+[startlam:nlams+startlam-1]);
end
totalpert=tert;
noise=nois;
injected=inj;
end