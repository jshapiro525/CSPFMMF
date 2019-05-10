
Zk=reshape(Z_actual(37,:),imdim,imdim);
Zkfm=reshape(Z(37,:),imdim,imdim);
myfig(Zk)
myfig(Zkfm)
myfig(reshape(delZ(37,:),imdim,imdim));

figure()
image(Zk(47:57,48:58),'CDataMapping','scaled')
figure()
image(Zkfm(47:57,48:58)*mean(mean(Zk(47:57,48:58)))/mean(mean(Zkfm(47:57,48:58))),'CDataMapping','scaled')
figure()
plot(Zk(52,48:58),'r')
hold on
plot(Zkfm(52,48:58),'b--')
legend('CSP Result','Forward Model')
xlabel('Pixel Location')
ylabel('Pixel Value')
hold off

figure()
image((Zk(47:57,48:58)),'CDataMapping','scaled')
figure()
image(Zkfm(47:57,48:58)*mean(mean(Zk(47:57,48:58)))/mean(mean(Zkfm(47:57,48:58))),'CDataMapping','scaled')
figure()
plot((Zk(47:57,53)),'r')
hold on
plot(Zkfm(47:57,53),'b--')
legend('CSP Result','Forward Model')
xlabel('Pixel Location')
ylabel('Pixel Value')
hold off