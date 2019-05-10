function [injtrunc, noisetrunc, ttrunc] = splituppic(injected, noise, totalpert)


        injbig=padarray(injected,[10 10],0,'both');
        injtrunc(:,:,:,1) = injbig(2:54,2:54,:);
        injtrunc(:,:,:,2) = injbig(2:54,35:87,:);
        injtrunc(:,:,:,3) = injbig(2:54,68:120,:);
        injtrunc(:,:,:,4) = injbig(35:87,2:54,:);
        injtrunc(:,:,:,5) = injbig(35:87,35:87,:);
        injtrunc(:,:,:,6) = injbig(35:87,68:120,:);
        injtrunc(:,:,:,7) = injbig(68:120,2:54,:);
        injtrunc(:,:,:,8) = injbig(68:120,35:87,:);
        injtrunc(:,:,:,9) = injbig(68:120,68:120,:);

        noisebig=padarray(noise,[10 10],0,'both');
        noisetrunc(:,:,:,1) = noisebig(2:54,2:54,:);
        noisetrunc(:,:,:,2) = noisebig(2:54,35:87,:);
        noisetrunc(:,:,:,3) = noisebig(2:54,68:120,:);
        noisetrunc(:,:,:,4) = noisebig(35:87,2:54,:);
        noisetrunc(:,:,:,5) = noisebig(35:87,35:87,:);
        noisetrunc(:,:,:,6) = noisebig(35:87,68:120,:);
        noisetrunc(:,:,:,7) = noisebig(68:120,2:54,:);
        noisetrunc(:,:,:,8) = noisebig(68:120,35:87,:);
        noisetrunc(:,:,:,9) = noisebig(68:120,68:120,:);

        tbig=padarray(totalpert,[10 10],0,'both');
        ttrunc(:,:,:,1) = tbig(2:54,2:54,:);
        ttrunc(:,:,:,2) = tbig(2:54,35:87,:);
        ttrunc(:,:,:,3) = tbig(2:54,68:120,:);
        ttrunc(:,:,:,4) = tbig(35:87,2:54,:);
        ttrunc(:,:,:,5) = tbig(35:87,35:87,:);
        ttrunc(:,:,:,6) = tbig(35:87,68:120,:);
        ttrunc(:,:,:,7) = tbig(68:120,2:54,:);
        ttrunc(:,:,:,8) = tbig(68:120,35:87,:);
        ttrunc(:,:,:,9) = tbig(68:120,68:120,:);
   
end