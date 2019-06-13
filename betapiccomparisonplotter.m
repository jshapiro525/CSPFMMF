
load('BPmfmodes')
imdim=75;


temp=sum(Z_true(37,:,1),3);
Zk=reshape(temp,imdim,imdim);

temp=sum(Z(37,:,1),3);
Zkfm=reshape(temp,imdim,imdim);

figure(1)
image(Zk(47:57,48:58),'CDataMapping','scaled')
figure(2)
image(Zkfm(47:57,48:58)*mean(mean(Zk(47:57,48:58)))/mean(mean(Zkfm(47:57,48:58))),'CDataMapping','scaled')
figure(3)
plot(Zk(52,48:58),'r')
hold on
plot(Zkfm(52,48:58),'b--')
legend('CSP Result','Forward Model')
xlabel('Pixel Location')
ylabel('Pixel Value')
hold off

figure(4)
image((Zk(47:57,48:58)),'CDataMapping','scaled')
figure(5)
image(Zkfm(47:57,48:58)*mean(mean(Zk(47:57,48:58)))/mean(mean(Zkfm(47:57,48:58))),'CDataMapping','scaled')
figure(6)
plot((Zk(47:57,53)),'r')
hold on
plot(Zkfm(47:57,53),'b--')
legend('CSP Result','Forward Model')
xlabel('Pixel Location')
ylabel('Pixel Value')
hold off