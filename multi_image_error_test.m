close all;
clear all;
clc;

factor = 32;

lams = 1.2:(1.8-1.2)/factor:1.8;
parangs = [1:31/factor:32];

for L = 2:length(lams)
    for P = 1:length(parangs)
        pos = (L-2)*length(parangs)+P;
        [n(pos), medianerror(pos), meanerror(pos)] = testFMimages(lams(1:L),parangs(1:P));
        %[n(pos), medianerror(pos), meanerror(pos)] = testFMimagesshort(lams(1:L),parangs(1:P));
        pos/(factor*(factor+1))*100
    end
end

figure()
plot(n,meanerror,'.')
xlabel('Number of Images')
ylabel('Mean Z Error')

figure()
plot(n,medianerror,'.')
xlabel('Number of Images')
ylabel('Median Z Error')

%linear fit
[medc1, inds] = polyfit(n,medianerror,1);
x = sort(n);
y = x*medc1(1) + medc1(2);
hold on
plot(x,y,'k--')
residualssquared1 = sum((y-medianerror(inds)).^2)

% 2nd order poly
[medc2, inds] = polyfit(n,medianerror,2);
[x,inds] = sort(n);
y2 = x.^2*medc2(1) + x*medc2(2) +medc2(3);
hold on
plot(x,y2,'r')



%3/2 fit
x = [[n.^(3/2)]' n' [n.^(1/2)]' ones(length(n),1)];
y = medianerror';
model = inv(x'*x)*x'*y;
hold on
y = model(1)*n.^(3/2) + model(2)*n+model(3)*n.^(1/2)+model(4);
plot(sort(n),y(inds),'b')


% % 3rd order poly
% [medc3, s] = polyfit(n,medianerror,3);
% [x,inds] = sort(n);
% y = x.^3*medc3(1) + x.^2*medc3(2) + x*medc3(3)+medc3(4);
% hold on
% plot(x,y,'r')
% residualssquared3 = sum((y-medianerror(inds)).^2)

legend('Data Points','Linear','2nd Order')%,'3rd Order')

print('imagenumbererror','-depsc')

