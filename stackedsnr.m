function [stack, SNRmap] = stackedsnr(images,parangs,lams)

N = length(lams); % Number of wavelengths in a datacube
dur = length(parangs); % How many datacubes you have

for i = 1:dur
    for j = 1:N
        pos = (i-1)*N+j; %Keeping track of which iteration you are across both loops
        aligned(:,:,pos) = imrotate(images(:,:,pos),parangs(i),'bicubic','crop');
    end
end

stack = sum(aligned,3);
squarestack = stack.^2;

myfigcbar(images(:,:,1))
myfigcbar(stack)
myfigcbar(squarestack)

end