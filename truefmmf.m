function [slicemfs] = truefmmf(Z, noise, injected1, parangs, lams, modes, currentpos)

imdim = size(noise,1);
stampsize = 21;
stampoff = floor(stampsize/2);
nlams = length(lams);
nmodes = size(Z,1);

slicemfs = zeros(imdim,imdim,modes,nlams);
parfor i = stampoff+1:imdim-stampoff-1
    for j = stampoff+1:imdim-stampoff-1
        
        injected = moveinjected2(injected1,currentpos,[i j],parangs,lams);
        
        [Z_est, delZ] = fmfinal(noise, injected, parangs, lams, modes);
        template = delZ;
        
        for lam = 1:nlams
            for mode = 1:modes
                k = nmodes-modes+mode;
                templatesquare = reshape(template(k,:,lam),imdim,imdim);
                datasquare = reshape(Z(k,:,lam),imdim,imdim);
                templatestamp = templatesquare(i-stampoff:i+stampoff,j-stampoff:j+stampoff);
                datastamp = datasquare(i-stampoff:i+stampoff,j-stampoff:j+stampoff);
                slicemfs(i,j,mode,lam) = sum(sum(datastamp.*templatestamp));
            end
        end
 
    end
end

end