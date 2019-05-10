function [alpha] = kstest(x1,cdf1,x2,cdf2)

n=length(cdf1);
m=length(cdf2);
D=max(abs(ecdf-pcdf));

c=D/(sqrt((n+m)/(n*m)));
alpha=exp(-2*c^2);

end
