clear;

[WIF,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZIF,GV,maxtau,PopS,CI] = LoadYemenData;


nd=WIF(GNZIF,(maxtau+1):end);
nd=length(nd(:));

load('Combo.mat');
beta=zeros(64,28);
DA=zeros(64,1);
K=zeros(2,64);
n=zeros(2,64);
dV=zeros(2,64);
KP=zeros(2,64);
r0=zeros(64,1);
KV=zeros(64,1);
w=zeros(64,1);
sigma_w=zeros(64,1);
DARv=zeros(64,1);
AIC=zeros(64,1);
k=zeros(64,1);
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
for ii=1:64
    load(['Fit-Vaccination-IncidenceperCapita' C(INC{ii}).N '-CalibratedDAR.mat']);
    [k(ii),beta(ii,:),tau,DB,DA(ii),K(:,ii),n(:,ii),KP(:,ii),KV(ii),dV(:,ii),r0(ii),~,w(ii),sigma_w(ii)]=RetParameterGA(par,XU,4);
    DARv(ii)=DAR;
    
    AIC(ii)=AICScore(k(ii),nd,RSSv);
end
K=K';
n=n';
KP=KP';
dV=dV';
AIC=AIC-min(AIC);
T=table(beta,DA,K,n,KP,KV,dV,r0,temp_0,DARv,w,sigma_w,AIC);
 
 writetable(T,'ModelParameter.csv','Delimiter',',');