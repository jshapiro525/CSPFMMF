
load('BPmfmodes')
imdim=75;



temp=sum(Z_true(37,:,1),3);

Zk=reshape(temp,imdim,imdim);


temp=sum(Z(37,:,1),3);
Zkfm=reshape(temp,imdim,imdim);

y=Zk(47:57,48:58);
x=Zkfm(47:57,48:58);

sum(sum(y.*x))/sum(sum(x.*x))
