close all

load('BPmfmodes')
imdim=75;


temp=sum(Z_true(37,:,1),3);
Zk=reshape(temp,imdim,imdim);

temp=sum(Z(37,:,1),3);
Zkfm=reshape(temp,imdim,imdim);

figure(1)
image(Zk(47:57,48:58),'CDataMapping','scaled')
axis square
figure(2)
image(Zkfm(47:57,48:58),'CDataMapping','scaled')
axis square

figure(3)
yyaxis left
plot(Zk(52,48:58),'r','LineWidth',3)
hold on
plot(Zkfm(52,48:58),'b--','LineWidth',3)
yyaxis right
plot(abs(Zkfm(52,48:58)-Zk(52,48:58)),'k-.','LineWidth',2)
legend('CSP Result','Forward Model','Difference','Fontsize',16)
xlabel('Pixel Location','Fontsize',16)
ylabel('Difference','Fontsize',16)
hold off
yyaxis left
ylabel('Pixel Value','Fontsize',16)
set(gca,'FontSize',16)

figure(4)
image((Zk(47:57,48:58)),'CDataMapping','scaled')
axis square
figure(5)
image(Zkfm(47:57,48:58),'CDataMapping','scaled')
axis square

figure(6)
yyaxis left
plot((Zk(47:57,53)),'r','LineWidth',2)
hold on
plot(Zkfm(47:57,53),'b--','LineWidth',2)
yyaxis right
plot(abs(Zkfm(47:57,53)-Zk(47:57,53)),'k-.','LineWidth',2)
legend('CSP Result','Forward Model','Difference','Fontsize',16)
xlabel('Pixel Location','Fontsize',16)
ylabel('Difference','Fontsize',16)
hold off
yyaxis left
ylabel('Pixel Value','Fontsize',16)
set(gca,'FontSize',16)