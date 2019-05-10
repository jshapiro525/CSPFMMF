function [result1, result2] = klipfmmfresults(scale)

filename=('BP_FMMF_PCA.fits');
BP=fitsread(filename,'image');
BPtrunc=BP(145-60:145+60,147-60:147+60);
BPtrunc=flipud(BPtrunc);
BProt=imrotate(BPtrunc,-120,'crop');
imdims=size(BProt);
result1=BProt(ceil(imdims(1)/2)-37:ceil(imdims(1)/2)+37,ceil(imdims(2)/2)-37:ceil(imdims(2)/2)+37);
if scale
    result1=result1-min(min(result1));
    result1=result1/max(max(result1))*100;
end

% myfigcbar(result1)

filename =('HD_FMMF_PCA.fits');
HD=fitsread(filename,'image');
HDtrunc=HD(145-60:145+60,146-60:146+60);
HDtrunc=flipud(HDtrunc);
HDrot=imrotate(HDtrunc,-80,'crop');
imdims=size(HDrot);
result2=HDrot(ceil(imdims(1)/2)-37:ceil(imdims(1)/2)+37,ceil(imdims(2)/2)-37:ceil(imdims(2)/2)+37);
if scale
    result2=result2-min(min(result2));
    result2=result2/max(max(result2))*100;
end
% myfigcbar(result2)

end