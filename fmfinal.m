function [Z, delZ] = fmfinal(noise, injected, parangs, lams, modes)

    [rep_parangs, rep_lams] = repparlam(parangs,lams);    
    
    for thislam =1:length(lams)
        % tperlam = totalpert(:,:,rep_lams==lams(thislam));
        nperlam = noise(:,:,rep_lams==lams(thislam));
        iperlam = injected(:,:,rep_lams==lams(thislam));
        [Z(:,:,thislam), Z_true(:,:,thislam), Zs(:,:,thislam), delZ(:,:,thislam), planetstack(:,:,thislam), fmstack(:,:,thislam)] = splitFM(nperlam, iperlam, modes, parangs,0);
    end


end
