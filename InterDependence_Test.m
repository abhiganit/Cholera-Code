%% Runs the mediation analysis for diesel
close all;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
% Dieslel
load('Fit-Vaccination-IncidenceperCapita-Diesel-CalibratedDAR.mat')
[k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
 
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XD=Xt((3.*maxtau+1):4*maxtau,:,:);
max(max(Xt(5,:,:)))
% Test for significance
RS=WI(GNZI,(maxtau+1):end)-Yt;
temp=squeeze(Xt(1,:,:));
XCv=[temp(:)];
for ii=2:28
    temp=squeeze(Xt(ii,:,:));
    XCv=[XCv temp(:)];
end
[SE,tStat,pValue] = LinRegressionStat(beta,RS(:),k,XCv,0);

% Conflict
load('Fit-Vaccination-IncidenceperCapita-Shellings-CalibratedDAR.mat')
[k,beta,tau,DB,DA,K,n,~,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
beta(5:8)
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
XC=Xt((1+maxtau):2.*maxtau,:,:);

% Test for significance
RS=WI(GNZI,(maxtau+1):end)-Yt;
temp=squeeze(Xt(1,:,:));
XCv=[temp(:)];
for ii=2:28
    temp=squeeze(Xt(ii,:,:));
    XCv=[XCv temp(:)];
end
[SE,tStat,pValue] = LinRegressionStat(beta,RS(:),k,XCv,0);


% First analysis

indx=[ 13:16 21:24 25:28];
% Test Dependence between Diesel and conflict
options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
[betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),-6.*ones(size(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
betap=10.^betap;
beta=zeros(1,28);
beta(indx)=betap;
[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

% Test for significance
RS=WI(GNZI,(maxtau+1):end)-Yt;
temp=squeeze(Xt(1,:,:));
XCv=[temp(:)];
for ii=2:28
    temp=squeeze(Xt(ii,:,:));
    XCv=[XCv temp(:)];
end
[SE,tStat,pValue] = LinRegressionStat(beta,RS(:),12,XCv,0);

fval=mean(residual.^2);

D=Xt(13,:,:);
D=D(:);
C=Xt(9,:,:);
C=C(:);
opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval1]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C)).^2)),log10([162 32]),[],[],[],[],[-32 -32],[3 3],[],opts);
bd=10.^bd;
bd2=bd;
[sum(bd(2).*beta(13:16))]
%save('DieselrepresentedthroughConflict.mat','bd','XC');
figure('units','normalized','outerposition',[0 0 1 1]);
scatter(C,D,40,'k','filled'); hold on
plot(linspace(min(C)*0.99,max(C)*1.01,10001),bd(1)+bd(2).*linspace(min(C)*0.99,max(C)*1.01,10001),'r','LineWidth',2);
xlabel('Targeted covariate','Fontsize',16);
ylabel('Diesel covariate','Fontsize',16);
box off;
[r,p]=corr(D,C);
text(1,max(ylim)*0.99,['r=' num2str(r) '(p=' num2str(p) ')'],'Fontsize',16);
text(1,max(ylim)*0.95,['D_t=' num2str(round(bd(1),2)) '+' num2str(round(bd(2),2)) 'C_t' ],'Fontsize',16);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',24);