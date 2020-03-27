function [] = graph_element_distro(matrix,Title,Y,X,L1)

[s1 s2] = size(matrix);
vect = reshape(matrix,s1*s2,1);

bins =linspace(log10(min(vect)+10e-15),log10(max(vect)),100);
bins = 10.^bins;
xc = histc(vect,bins);

figure()
semilogx(bins,xc,'r','LineWidth',2.5)

hold on
h = line([median(vect) median(vect)],[0 max(xc)],'LineWidth',2.5);
L = line([mean(vect) mean(vect)],[0 max(xc)],'color','K','LineStyle','--' ,'LineWidth',2.5);
title(Title)
ylabel(Y,'FontSize',18)
xlabel(X,'FontSize',18)
legend(L1,horzcat(['Median = ' num2str(median(vect)) '%']),horzcat(['Mean = ' num2str(mean(vect)) '%']),'FontSize',18)
axis([10^-4 10^4 0 6*10^5])

end