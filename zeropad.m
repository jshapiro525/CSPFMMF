function [new] = zeropad(x,n)

nrows=size(x,1);
ncols=size(x,2);
x=[zeros(nrows,n) x zeros(nrows,n)];
new=[zeros(n,2*n+ncols); x; zeros(n,2*n+ncols)];

end