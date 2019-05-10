close all
yellow = 5+rand(2,2)*2-1;
green = 3+rand(2,2)*2-1;
blue = rand(2,2)*2-1;

pos2=[1 1 2 2];
pos1=[1 2 1 2];

spaces='ABCD';

for i=1:3
    
    if i==1
        myfig(yellow)
    elseif i==2
        myfig(green)
    else
        myfig(blue)
    end
    for j=1:4
        text(pos1(j),pos2(j),spaces(j),'FontSize',14,'HorizontalAlignment','center')
    end
    axis off
    caxis([0 max(max(yellow))])
end

for i=1:4
    