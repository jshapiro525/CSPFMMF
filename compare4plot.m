load('Final4res.mat')

BPCSP=BPCSP./(max(max(BPCSP)))*100;
HDCSP=HDCSP./(max(max(HDCSP)))*100;
BPKLIP=BPKLIP./(max(max(BPKLIP)))*100;
HDKLIP=HDKLIP./(max(max(HDKLIP)))*100;

figure(1)
imagesc(BPCSP)
axis square
axis off
colorbar
caxis([-90 100])

figure(2)
imagesc(HDCSP)
axis square
axis off
colorbar
caxis([-90 100])

figure(3)
imagesc(BPKLIP)
axis square
axis off
colorbar
caxis([-90 100])

figure(4)
imagesc(HDKLIP)
axis square
axis off
colorbar
caxis([-90 100])