function [] = log_graph_element_distro(matrix,Title,Y,X,L1)

[s1, s2] = size(matrix);
vect = reshape(matrix,s1*s2,1);

bins =linspace(log10(min(vect)+1e-16),log10(max(vect)),100);
bins = 10.^bins;
xc = histc(vect,bins);

figure()
semilogx(bins,xc,'r')
yl=ylim;
hold on
h = line([median(vect) median(vect)],[0 max(yl)]);
L = line([mean(vect) mean(vect)],[0 max(yl)],'color','k','LineStyle','--');
title(Title)
ylabel(X)
xlabel(Y)
legend(L1,horzcat(['Median = ' num2str(median(vect)) '%']),horzcat(['Mean = ' num2str(mean(vect)) '%']))

end