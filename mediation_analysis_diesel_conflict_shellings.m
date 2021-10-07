%% Runs the mediation analysis for conflict & shellings via diesel
close all;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
% Dieslel
load('Fit-Vaccination-IncidenceperCapita-Diesel-Rain-CalibratedDAR.mat')
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
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Rain-CalibratedDAR.mat')
[k,beta,tau,DB,DA,K,n,~,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);

[Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);


% Test for significance
RS=WI(GNZI,(maxtau+1):end)-Yt;
temp=squeeze(Xt(1,:,:));
XCv=[temp(:)];
for ii=2:28
    temp=squeeze(Xt(ii,:,:));
    XCv=[XCv temp(:)];
end
[SE,tStat,pValue] = LinRegressionStat(beta,RS(:),k,XCv,0);
clc;
beta(5:12)
pValue(5:12)


% First analysis


indx=[5:8 9:12 13:16 21:24 25:28];
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
beta(5:12)
pValue(5:12)
fval=mean(residual.^2);

XC=Xt((1+maxtau):2.*maxtau,:,:);
XS=Xt((1+2.*maxtau):3.*maxtau,:,:);

D=Xt(13,:,:);
D=D(:);


C=XC(1,:,:);
C=C(:);
C2=XC(2,:,:);
C2=C2(:);

S=XS(1,:,:);
S=S(:);
opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval1]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C+10.^x(3).*S+10.^x(4).*C2)).^2)),log10([162 31 1 1 ]),[],[],[],[],[-32 -32 -32 -32],[3 3 3 3],[],opts);
bd=10.^bd;
bd2=bd;
[sum(bd(3).*beta(13:16))]
[sum(bd(2).*beta(5:8))]
%save('DieselrepresentedthroughConflictShellings.mat','bd','XC','XS');

% Second analysis
load('Fit-Vaccination-IncidenceperCapita-Diesel-Rain-CalibratedDAR.mat')
[~,beta,tau,DB,DA,~,~,KP,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
indx=[5:8 9:12 13:16 21:24 25:28];
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
beta(5:12)
pValue(5:12)

XC=Xt((1+maxtau):2.*maxtau,:,:);
XS=Xt((1+2.*maxtau):3.*maxtau,:,:);

fval=mean(residual.^2);
D=Xt(13,:,:);
D=D(:);

C=XC(1,:,:);
C=C(:);
C2=XC(2,:,:);
C2=C2(:);

S=XS(1,:,:);
S=S(:);
opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval2]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C+10.^x(3).*S+10.^x(4).*C2)).^2)),log10([69 18 0.01 12 ]),[],[],[],[],[-32 -32 -32 -32],[3 3 3 3],[],opts);
bd=10.^bd;
bd2=bd;
[sum(bd(3).*beta(13:16))]
[sum(bd(2).*beta(5:8))]
save('DieselrepresentedthroughConflictShellings.mat','bd','XC','XS');
