function [injected]=injectplanets(base,i,j,imdim,maxsigstr)

x1 = 1:1:imdim; x2 = 1:1:imdim;
[X1,X2] = meshgrid(x1,x2);

injected=base;

for m=1:length(i)
    mu=[i(m) j(m)];
    Sigma = [2 0; 0 2];
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    F = F./max(max(F))*maxsigstr;
    injected=injected + F;
end

end
