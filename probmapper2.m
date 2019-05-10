function [p1,esnr1,lambda] = probmapper2(val, vect)

skew=skewness(vect);
sigma=std(vect);
kurt=kurtosis(vect);
mu=mean(vect);


var=sigma^2;
place=0;
nunits=1e4;
distrib=[linspace(-2*val,val,nunits) linspace(val,2*val,nunits/10)];

disp('Doing PearsPDF')

for i=distrib
    place=place+1;  
    [prob(place),type,~,lambda] = pearspdf(i,mu,sigma,skew,kurt);
end

pcdf=cumtrapz(prob);
normfactor=max(pcdf);
pcdf=pcdf/normfactor;
prob=prob/normfactor;


disp('Starting CDF Integration')
tempsum=0;
spots=[nunits/4 nunits/2 nunits/4*3];

for i=1:nunits-1
    tempsum=tempsum+(prob(i)+prob(i+1))/2;
    if ismember(i,spots)
        i/(nunits-1)
    end
end

p1=tempsum;
esnr1=abs(double(norminv(1-tempsum)));

figure
semilogy(distrib,prob)
hold on
plot(distrib, normpdf(distrib,mu,sigma))
plot(val,prob(nunits),'x','MarkerSize',20)
plot(max(vect),pearspdf(max(vect),mu,sigma,skew,kurt)/normfactor,'o','MarkerSize',12)
legend('4-Norm Distribution','Normal Distribution','Value','Maximum Empirical Data','Location','southeast')
axis([-2*val 2*val 10^-20 1])

[M,I]=min(abs(normpdf(distrib,mu,sigma)-prob(nunits)));
pos=distrib(I);
% esnr2=abs((pos-mu)/sigma);

% p1=pcdf(nunits);
% esnr1=norminv(p1);

vectedsort=sort(vect);
ecdf=[];
for t=1:length(vectedsort)
    ecdf(t)=t/length(vectedsort);
end
figure()
plot(vectedsort,ecdf,'.',distrib,pcdf,'-');
axis([-val/2 val+10 0 1])
legend('Empirical CDF','Pearson Distribution Fit CDF','Location','southeast')

end



