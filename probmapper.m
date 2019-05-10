function [cdfval, esnr, lambda ] = probmapper(val, vect)

skew=skewness(vect);
sigma=std(vect);
kurt=kurtosis(vect);
mu=mean(vect);
lambda=mu;

var=sigma^2;

Beta1=skew^2;
Beta2=kurt;

b0=(4*Beta2-3*Beta1)/(10*Beta2-12*Beta1-18)*var;
b1=sqrt(var)*sqrt(Beta1)*(Beta2+3)/(10*Beta2-12*Beta1-18);
b2=(2*Beta2-3*Beta1-6)/(10*Beta2-12*Beta1-18);

%y = (x+b1/(2*b2));

if b1^2-4*b2*b0<0
    alpha=sqrt(4*b2*b0-b1^2)/(2*b2);
    m=.5/b2;
    v=-(2*b2*b1+b1)/(2*b2^2*alpha);    
    
    lognorm = 2*log(abs(gammac(m+(v/2)*i)/gammac(m)))-log(alpha*beta(m-.5,.5));
    
    
    p=@(x) (  exp(lognorm + log((1+((x-lambda)./alpha).^2).^-m.*exp(-v.*atan((x-lambda)./alpha))))  );
   
    
else
        disp('Not Type IV distribution. Quitting.')
        return
end


nunits=1e3;
distrib=[linspace(-2*val,val,nunits) linspace(val,2*val,nunits/10)];

prob = p(distrib);

figure
semilogy(distrib,prob)
hold on
plot(distrib, normpdf(distrib,mu,sigma))
plot(val,prob(nunits),'x','MarkerSize',20)
plot(max(vect),p(max(vect)),'o','MarkerSize',12)
legend('4-Norm Distribution','Normal Distribution','Value','Maximum Empirical Data','Location','southeast')
axis([-2*val 2*val 10^-20 1])

for j=1:length(distrib)
    pcdf(j)=integral(p,-inf,distrib(j));
end

vectedsort=sort(vect);
ecdf=[];

for t=1:length(vectedsort)
    ecdf(t)=t/length(vectedsort);
end

figure()
plot(vectedsort,ecdf,'.',distrib,pcdf,'-');
axis([-val/2 val+10 0 1])
legend('Empirical CDF','Pearson Distribution Fit CDF','Location','southeast')





pval=p(distrib(nunits));
cdfval=pcdf(nunits);
esnr=norminv(pcdf(nunits));


end



