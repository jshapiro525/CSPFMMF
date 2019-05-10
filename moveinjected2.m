function [newA] = moveinjected2(A,currentpos,newpos,parangs,lams)

n=length(A(1,1,:));
[rep_parangs, rep_lams] = repparlam(parangs,lams);
y0=currentpos(1);
x0=currentpos(2);
imdim = size(A(:,:,1),1);
y2 = newpos(1);
x2 = newpos(2);


for i=1:n
    newderotA = imrotate(A(:,:,i),rep_parangs(i)-parangs(1),'bicubic','crop');
    newposA = imtranslate(newderotA,[x2-x0 y2-y0]);
    newA(:,:,i) = imrotate(newposA,-1*(rep_parangs(i)-parangs(1)),'bicubic','crop');
end

end
