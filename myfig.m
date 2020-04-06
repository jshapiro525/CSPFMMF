function [] = myfig(x)

figure()
image(x,'cdatamapping','scaled')
colorbar
    hAxes = gca;
    axis(hAxes,'square')
% %    hAxes.XRuler.Axle.LineStyle = 'none';  
   axis off

end
