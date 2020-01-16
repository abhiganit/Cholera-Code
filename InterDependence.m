%% Runs the mediation analysis for diesel
close all;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
load('Fit-Vaccination-IncidenceperCapita-Conflict-CalibratedDAR.mat')
[~,beta,tau,DB,DA,K2,n2,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
betaN=beta;
%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K2,n2,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XC=Xt((1+maxtau):2.*maxtau,:,:);


load('Fit-Vaccination-IncidenceperCapita-Shellings-CalibratedDAR.mat')
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
K2(2)=K(2);
n2(2)=n(2);
%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,Xt]=  LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XS=Xt((1+2.*maxtau):3.*maxtau,:,:);

load('Fit-Vaccination-IncidenceperCapita-Diesel-CalibratedDAR.mat')
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XD=Xt((3.*maxtau+1):4*maxtau,:,:);

indx=[5:8];
% Test Dependence between Diesel and conflict
options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
[betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K2,n2,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),-6.*ones(size(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
betap=10.^betap;
fval=mean(residual.^2);
temp=squeeze(XC(1,:,:));
XCv=[temp(:)];
temp=squeeze(XC(2,:,:));
XCv=[XCv temp(:)];
temp=squeeze(XC(3,:,:));
XCv=[XCv temp(:)];
temp=squeeze(XC(4,:,:));
XCv=[XCv temp(:)];

temp=squeeze(XS(1,:,:));
XSv=[temp(:)];
temp=squeeze(XS(2,:,:));
XSv=[XSv temp(:)];
temp=squeeze(XS(3,:,:));
XSv=[XSv temp(:)];
temp=squeeze(XS(4,:,:));
XSv=[XSv temp(:)];

test=RSSv;
[SE,tStat,pValue] = LinRegressionStat(betap,residual,4,XCv,0);


D=XD(1,:,:);
D=D(:);
C=XC(1,:,:);
C=C(:);
S=XS(1,:,:);
S=S(:);
opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval1]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C)).^2)),log10([162 32]),[],[],[],[],[-32 -32],[3 3],[],opts);
bd=10.^bd;
bd2=bd;
[sum(bd(2).*beta(13:16))]
save('DieselrepresentedthroughConflict.mat','bd','XC');
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(C,D,40,'k','filled'); hold on
plot(linspace(min(C)*0.99,max(C)*1.01,10001),bd(1)+bd(2).*linspace(min(C)*0.99,max(C)*1.01,10001),'r','LineWidth',2);
xlabel('Conflict covariate','Fontsize',16);
ylabel('Diesel covariate','Fontsize',16);
box off;
[r,p]=corr(D,C);
text(1,max(ylim)*0.99,['r=' num2str(r) '(p=' num2str(p) ')'],'Fontsize',16);
text(1,max(ylim)*0.95,['D_t=' num2str(round(bd(1),2)) '+' num2str(round(bd(2),2)) 'C_t' ],'Fontsize',16);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',24);
[SE,tStat,pValue] = LinRegressionStat(bd,(D-(bd(1)+bd(2).*S)),2,[ones(1,length(S)); S']',0);

test=0;

indx=[9:12];
% Test Dependence between Diesel and conflict
options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
[betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K2,n2,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),-6.*ones(size(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
betap=10.^betap;
fval=mean(residual.^2);
temp=squeeze(XC(1,:,:));
XCv=[temp(:)];
temp=squeeze(XC(2,:,:));
XCv=[XCv temp(:)];
temp=squeeze(XC(3,:,:));
XCv=[XCv temp(:)];
temp=squeeze(XC(4,:,:));
XCv=[XCv temp(:)];

temp=squeeze(XS(1,:,:));
XSv=[temp(:)];
temp=squeeze(XS(2,:,:));
XSv=[XSv temp(:)];
temp=squeeze(XS(3,:,:));
XSv=[XSv temp(:)];
temp=squeeze(XS(4,:,:));
XSv=[XSv temp(:)];

test=RSSv;
[SE,tStat,pValue] = LinRegressionStat(betap,residual,4,XSv,0);


D=XD(1,:,:);
D=D(:);
C=XC(1,:,:);
C=C(:);
S=XS(1,:,:);
S=S(:);
opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval2]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*S)).^2)),log10([162 32]),[],[],[],[],[-32 -32],[3 3],[],opts);
bd=10.^bd;
[sum(bd(2).*beta(13:16))]
bd2=[(bd(1)+bd2(1))./2 bd2(2) bd(2)];
save('DieselrepresentedthroughShellings.mat','bd','XC');
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(S,D,40,'k','filled'); hold on
plot(linspace(min(S)*0.99,max(S)*1.01,10001),bd(1)+bd(2).*linspace(min(S)*0.99,max(S)*1.01,10001),'r','LineWidth',2);
xlabel('Shellings covariate','Fontsize',16);
ylabel('Diesel covariate','Fontsize',16);
box off;
[r,p]=corr(D,S);
text(1,max(ylim)*0.99,['r=' num2str(r) '(p=' num2str(p) ')'],'Fontsize',16);
text(1,max(ylim)*0.95,['D_t=' num2str(round(bd(1),2)) '+' num2str(round(bd(2),2)) 'S_t' ],'Fontsize',16);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',24);

[SE,tStat,pValue] = LinRegressionStat(bd,(D-(bd(1)+bd(2).*S)),2,[ones(1,length(S)); S']',0);

test=0;

[bd,fval3]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C+10.^x(3).*S)).^2)),[2.1 1.2 0.2],[],[],[],[],[-32 -32 -32],[3 3 3],[],opts);
bd=10.^bd;
save('DieselrepresentedthroughConflictShellings.mat','bd','XC','XS');


%% District

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData;
load('Fit-Vaccination-IncidenceperCapita-Conflict-CalibratedDAR.mat')
[~,beta,tau,DB,DA,K2,n2,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
betaN=beta;
n2=n2-TruncV;
n2(n2<1)=1;

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K2,n2,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XC=Xt((1+maxtau):2.*maxtau,:,:);


load('Fit-Vaccination-IncidenceperCapita-Shellings-CalibratedDAR.mat')
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
n=n-TruncV;
n(n<1)=1;

K2(2)=K(2);
n2(2)=n(2);
%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,Xt]=  LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XS=Xt((1+2.*maxtau):3.*maxtau,:,:);

save('DieselrepresentedthroughConflictShellings_District.mat','bd','XC','XS');

