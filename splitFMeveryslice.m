function [Z_true, Phi_true, Z, Zs, delZ] = splitFMeveryslice(totalpert, noise, injected,parangs, lams, modes, figs)

    [rep_parangs, rep_lams] = repparlam(parangs,lams);    
    
    for thislam =1:length(lams)
        tperlam = totalpert(:,:,rep_lams==lams(thislam));
        iperlam = injected(:,:,rep_lams==lams(thislam));
        nperlam = noise(:,:,rep_lams==lams(thislam));
        [Z_true(:,:,thislam), Phi_true(:,:,thislam), Z(:,:,thislam), Zs(:,:,thislam), delZ(:,:,thislam)] = splitFM(nperlam, iperlam, modes, parangs,figs);
        fprintf(horzcat(['FM is ' num2str(thislam/length(lams)*100) ' percent done\n']))               
    end


end
