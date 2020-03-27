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
plot(n,abs(medianerror-1),'.')
xlabel('Number of Images')
ylabel('Normalized Photometric Error')



print('imagenumbererror','-depsc')

