load('Fit-Vaccination-PercentData=80.mat');
%clear;
%clc;
% load('Fit-Vaccination-PercentData=80-Rainfall.mat')
AIC=zeros(length(par(:,1)),1);
BIC=zeros(length(par(:,1)),1);
PDS=0.8;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,PDS);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
nd=WI(GNZI(GTF),1:NW);
nd=length(nd(:));
for ii=1:length(par(:,1))
    XUt=XU;
    [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,mln,a,KV,dV]=RetParameterPS(par(ii,:),XUt,[2;2]);
    XUt(beta<10^(-4))=0;
    beta(XU==1)
    [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,mln,a,KV,dV]=RetParameterPS(par(ii,:),XUt,[2;2]);
    k
    
    AIC(ii)=AICScore(k,nd,RSSv(ii)*nd);
    BIC(ii)=BICScore(k,nd,RSSv(ii)*nd);
end
% RSSv
% CVE
% dAIC=AIC-min(AIC)


