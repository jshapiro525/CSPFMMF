function [slicemfs] = segmentedtruefmmf(Z, noise, injected1, totalpert, parangs, lams, modes, currentpos)

imdim = 53;
stampsize = 21;
stampoff = floor(stampsize/2);
nlams = length(lams);
nmodes = size(Z,1);
currentpos=currentpos+9;

slicemfstrunc = zeros(imdim,imdim,modes,nlams,9);
parfor i = stampoff+1:imdim-stampoff-1
    for j = stampoff+1:imdim-stampoff-1
        
        injected = moveinjected2(injected1,currentpos,[i j],parangs,lams);
        
        [injtrunc, noisetrunc, ttrunc] = splituppic(injected, noise, totalpert);

        for K=1:9
            
            [Z_est, delZ] = fmfinal(noisetrunc(:,:,:,K), injtrunc(:,:,:,K), parangs, lams, modes);
            template = delZ;

            for lam = 1:nlams
                for mode = 1:modes
                    k = nmodes-modes+mode;
                    templatesquare = reshape(template(k,:,lam),imdim,imdim);
                    datasquare = reshape(Z(k,:,lam,K),imdim,imdim);
                    templatestamp = templatesquare(i-stampoff:i+stampoff,j-stampoff:j+stampoff);
                    datastamp = datasquare(i-stampoff:i+stampoff,j-stampoff:j+stampoff);
                    slicemfstrunc(i,j,mode,lam,K) = sum(sum(datastamp.*templatestamp));
                end
            end 
        end
    end
end

slicemfs1 = cat(2,slicemfstrunc(11:43,11:43,:,:,1),slicemfstrunc(11:43,11:43,:,:,2),slicemfstrunc(11:43,11:43,:,:,3));
slicemfs2 = cat(2,slicemfstrunc(11:43,11:43,:,:,4),slicemfstrunc(11:43,11:43,:,:,5),slicemfstrunc(11:43,11:43,:,:,6));
slicemfs3 = cat(2,slicemfstrunc(11:43,11:43,:,:,7),slicemfstrunc(11:43,11:43,:,:,8),slicemfstrunc(11:43,11:43,:,:,9));
slicemfs = cat(1,slicemfs2,slicemfs2,slicemfs3);

end