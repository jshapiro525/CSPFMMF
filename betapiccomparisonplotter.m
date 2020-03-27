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
axis off

figure(2)
image(Zkfm(47:57,48:58),'CDataMapping','scaled')
axis square
axis off

figure(3)
yyaxis left
plot(Zk(52,48:58),'r','LineWidth',3)
hold on
plot(Zkfm(52,48:58),'b--','LineWidth',3)
yyaxis right
plot(abs(Zkfm(52,48:58)-Zk(52,48:58)),'k-.','LineWidth',3)
legend('CSP Result','Forward Model','Difference','Fontsize',20)
xlabel('Pixel Location','Fontsize',20)
ylabel('Difference','Fontsize',20)
hold off
yyaxis left
ylabel('Pixel Value','Fontsize',20)
set(gca,'FontSize',20)

figure(4)
image((Zk(47:57,48:58)),'CDataMapping','scaled')
axis square
axis off

figure(5)
image(Zkfm(47:57,48:58),'CDataMapping','scaled')
axis square
axis off

figure(6)
yyaxis left
plot((Zk(47:57,53)),'r','LineWidth',3)
hold on
plot(Zkfm(47:57,53),'b--','LineWidth',3)
yyaxis right
plot(abs(Zkfm(47:57,53)-Zk(47:57,53)),'k-.','LineWidth',3)
legend('CSP Result','Forward Model','Difference','Fontsize',20)
xlabel('Pixel Location','Fontsize',20)
ylabel('Difference','Fontsize',20)
hold off
yyaxis left
ylabel('Pixel Value','Fontsize',20)
set(gca,'FontSize',20)
