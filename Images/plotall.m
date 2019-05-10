function [] = plotall(images)
n=length(images(:,1,1));
LAM=1;
for i=1:n
%    subplot(ceil(n/5),5,i)
    figure()
    image(reshape(images(i,:,LAM),75,75),'cdatamapping','scaled')

    hAxes = gca;
    axis(hAxes,'square')
%    hAxes.XRuler.Axle.LineStyle = 'none';  
    axis off
    title(num2str(i))
 %   fname=horzcat(['Zmodes',num2str(i)]);
%    print(fname,'-depsc')
end
end
