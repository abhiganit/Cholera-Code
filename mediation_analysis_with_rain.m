%% Runs the mediation analysis for diesel
clear;
close all;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

ModelBase={'Diesel','Conflict','Conflict-Diesel','Diesel','Shellings','Shellings-Diesel','Wheat','Conflict','Conflict-Wheat','Wheat','Shellings','Shellings-Wheat'};

CovariateConflict={'Conflict','Conflict','Conflict','Shellings','Shellings','Shellings','Conflict','Conflict','Conflict','Shellings','Shellings','Shellings'};

CovariatePrice={'Diesel','Diesel','Diesel','Diesel','Diesel','Diesel','Wheat','Wheat','Wheat','Wheat','Wheat','Wheat'};

indxv=[5:8 13:16 21:28; 5:8 13:16 21:28; 5:8 13:16 21:28; 9:12 13:16 21:28; 9:12 13:16 21:28; 9:12 13:16 21:28; 5:8 17:20 21:28; 5:8 17:20 21:28; 5:8 17:20 21:28; 9:12 17:20 21:28; 9:12 17:20 21:28; 9:12 17:20 21:28;];

Boot_betaC=cell(length(ModelBase),1);
Boot_betaC2=cell(length(ModelBase),1);
pValue_beta_Conflict=zeros(length(ModelBase),4);
betaC=zeros(length(ModelBase),4);
betaC2=zeros(length(ModelBase),4);
pValue_beta_Mediation=zeros(length(ModelBase),4);
pValue_Corr=zeros(length(ModelBase),2);
betaCorr=zeros(length(ModelBase),2);
for mm=1:length(ModelBase)

    % Load the baseline model
    load(['Fit-Vaccination-IncidenceperCapita-' ModelBase{mm} '-Rain-CalibratedDAR.mat'],'par','XU','CF','maxtau','DAR','RF')
    [k,~,tau,~,~,~,~,~,KV,dV,r0,~,w]=RetParameterGA(par,XU,CF,maxtau);
    % Load Coefficient for conflict to see if there is a significinat
    % relationship
    load(['Fit-Vaccination-IncidenceperCapita-' CovariateConflict{mm} '-Rain-CalibratedDAR.mat'],'par','XU','CF','maxtau','RF')
    [~,beta,~,DB,DA,K,n,~,~,~,~,~,~]=RetParameterGA(par,XU,CF,maxtau);
    % Load the parameters for the prices    
    load(['Fit-Vaccination-IncidenceperCapita-' CovariatePrice{mm} '-Rain-CalibratedDAR.mat'],'par','XU','CF','maxtau','RF')
    [~,~,~,~,~,~,~,KP,~,~,~,~,~]=RetParameterGA(par,XU,CF,maxtau);
    if(strcmp(CovariateConflict{mm},ModelBase{mm}))
        [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
    else
        indx=indxv(mm,[1:4 9:16]);
        options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,zeros(1,28),tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),-6.*ones(size(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,28);
        beta(indx)=betap;
        [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
    end
    % Test for significance
    RS=WI(GNZI,(maxtau+1):end)-Yt;
    temp=squeeze(Xt(1,:,:));
    XCv=[temp(:)];
    for ii=2:28
        temp=squeeze(Xt(ii,:,:));
        XCv=[XCv temp(:)];
    end
    [SE,tStat,pValue] = LinRegressionStat(beta,RS(:),length(indx),XCv,0);

    betaCIC=beta;
    pValue_beta_Conflict(mm,:)=pValue(indxv(mm,1:4));
    betaC(mm,:)=[beta(indxv(mm,1:4))];


    % First analysis
    % Use the estimates from the conflict analysis but the threshold estimates
    % for the diese price analysis

    indx=indxv(mm,:);
    % Test Dependence between Diesel and conflict
    options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
    [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),-6.*ones(size(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
    betap=10.^betap;
    beta=zeros(1,28);
    beta(indx)=betap;

    betaMIC=beta;

    betaC2(mm,:)=[beta(indxv(mm,1:4))];
    [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

    % Test for significance
    RS=WI(GNZI,(maxtau+1):end)-Yt;
    temp=squeeze(Xt(1,:,:));
    XCv=[temp(:)];
    for ii=2:28
        temp=squeeze(Xt(ii,:,:));
        XCv=[XCv temp(:)];
    end
    [SE,tStat,pValue] = LinRegressionStat(beta,RS(:),length(indx),XCv,0);
    
    pValue_beta_Mediation(mm,:)=pValue(indxv(mm,1:4));
    D=Xt(indxv(mm,5),:,:);
    D=D(:);
    C=Xt(indxv(mm,1),:,:);
    C=C(:);
    opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
    [bd,fval1]=fmincon(@(x)(mean((D-(10.^x(1)+10.^x(2).*C)).^2)),log10([162 32]),[],[],[],[],[-1 -1],[3 3],[],opts);

    RS=D-(10.^bd(1)+10.^bd(2).*C);
    temp=ones(size(C));
    XCv=[temp(:) C];    
    [SE,tStat,pValue_Corr(mm,:)] = LinRegressionStat(10.^bd,RS(:),2,XCv,0);
    betaCorr(mm,:)=10.^bd;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Bootstrap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    betaCv=zeros(1000,4);
    betaC2v=zeros(1000,4);
    
        options = optimoptions('lsqnonlin','FunctionTolerance',10^(-10),'StepTolerance',10^(-10));
    parfor ii=1:1000
        % conflict only
        GNZIt=GNZI(randi(length(GNZI),size(GNZI)));
        indx=indxv(mm,[1:4 9:16]);
        % Test Dependence between Diesel and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaCIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,CF,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),RF,r0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaCIC(betaCIC~=0)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,28);
        beta(indx)=betap;
    
        betaCv(ii,:)=beta(indxv(mm,1:4)); 
        
        
        indx=indxv(mm,:);
        % Test Dependence between Diesel and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaMIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,CF,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),RF,r0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaMIC(betaMIC~=0)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,28);
        beta(indx)=betap;
    
        betaC2v(ii,:)=beta(indxv(mm,1:4)); 
    end
    Boot_betaC{mm}=betaCv;
    Boot_betaC2{mm}=betaC2v;
    save('Mediation_Analysis_Rainfall.mat');
end