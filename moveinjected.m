function [newA] = moveinjected(A,currentpos,newpos,parangs,lams)

n=length(A(1,1,:));
[rep_parangs, rep_lams] = repparlam(parangs,lams);
y0=currentpos(1);
x0=currentpos(2);
imdim = size(A(:,:,1),1);
xc = ceil(imdim/2);
yc = ceil(imdim/2);
r0 = sqrt((x0-xc)^2+(y0-yc)^2);
th0 = atan2(yc-y0,x0-xc)*180/pi;
y2 = newpos(1);
x2 = newpos(2);
r2 = sqrt((x2-xc)^2+(y2-yc)^2);
th2 = atan2(yc-y2,x2-xc)*180/pi;

thr=(th2-th0);

close all

for i=1:n
    newderotA = imrotate(A(:,:,i),rep_parangs(i)-parangs(1),'bicubic','crop');
    newposA = imtranslate(newderotA,[x2-x0 y2-y0]);
    newA(:,:,i) = imrotate(newposA,-1*(rep_parangs(i)-parangs(1)),'bicubic','crop');
end

end
