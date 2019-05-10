function [result1, result2] = klipfmmfresultssaved(scale)


load('KLIPdata.mat')

if scale
    result1=result1-min(min(result1));
    result1=result1/max(max(result1))*100;

    result2=result2-min(min(result2));
    result2=result2/max(max(result2))*100;
end


end