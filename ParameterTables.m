clear;

[WIF,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZIF,GV,maxtau,PopS,CI] = LoadYemenData;


nd=WIF(GNZIF,(maxtau+1):end);
nd=length(nd(:));

load('Combo.mat');
beta=zeros(length(INC),32);
DA=zeros(length(INC),1);
K=zeros(2,length(INC));
n=zeros(2,length(INC));
dV=zeros(2,length(INC));
KP=zeros(2,length(INC));
r0=zeros(length(INC),1);
temp_0=zeros(length(INC),1);
KV=zeros(length(INC),1);
w=zeros(length(INC),1);
sigma_w=zeros(length(INC),1);
DARv=zeros(length(INC),1);
AIC=zeros(length(INC),1);
k=zeros(length(INC),1);
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});
for ii=1:length(INC)
    load(['Fit-Vaccination-IncidenceperCapita' C(INC{ii}).N '-CalibratedDAR.mat']);
    [k(ii),beta(ii,:),tau,DB,DA(ii),K(:,ii),n(:,ii),KP(:,ii),KV(ii),dV(:,ii),r0(ii),temp_0(ii),~,w(ii),sigma_w(ii)]=RetParameter_Table(par,XU,4);
    DARv(ii)=DAR;
    
    AIC(ii)=aicbic(-RSSv,k(ii));
end
K=K';
n=n';
KP=KP';
dV=dV';
AIC=AIC-min(AIC);
T=table(beta,DA,K,n,KP,KV,dV,r0,temp_0,DARv,w,sigma_w,AIC);
 
 writetable(T,'ModelParameter.csv','Delimiter',',');