%% Runs the mediation analysis for diesel
clear;
close all;
clc;
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

load('Combo.mat');
CName={'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'};

% Models without Diesel but with Conflict or Shellings
ModelUD={};
ModelU={};
ModelConflictDiesel={};
ModelConflict={};
count=0;
for ii=1:length(INC)
    test=INC{ii};
    if(ismember(5,test) & (ismember(2,test)| ismember(3,test)))
        count=count+1;
        ModelUD{count}=[test];
        lt=[];
        ltc=[];
        tt=[];
        for jj=1:length(test)
            lt=[lt 4.*(test(jj)-1)+[1:4]];
            if(test(jj)~=5)
                tt=[tt test(jj)];
                ltc=[ltc 4.*(test(jj)-1)+[1:4]];
            end
        end
        
        ModelU{count}=[tt];
        lt=[lt 21:32];
        ltc=[ltc 21:32];
        ModelConflictDiesel{count}=lt;
        ModelConflict{count}=ltc;
    end
end

Boot_betaC=cell(length(ModelConflict),1);
Boot_betaC2=cell(length(ModelConflict),1);
pValue_beta_Conflict=zeros(length(ModelConflict),8);
betaC=zeros(length(ModelConflict),8);
betaC2=zeros(length(ModelConflict),8);
pValue_beta_Mediation=zeros(length(ModelConflict),8);
pValue_Corr=ones(length(ModelConflict),4);
betaCorr=zeros(length(ModelConflict),4);
for mm=1:length(ModelConflict)
       
    % Load the baseline model with Diesel and Conflict
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelUD{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [kD,betaM,~,~,~,~,~,KPd,~,~,~,~,~,~]=RetParameterGA(par,XU,maxtau);
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelU{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,sigma_w]=RetParameterGA(par,XU,maxtau);
    KP(2)=KPd(2);
    [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
    % Test for significance
    RS=WI(GNZI,(maxtau+1):end)-Yt;
    temp=squeeze(Xt(1,:,:));
    XCv=[temp(:)];
    for ii=2:32
        temp=squeeze(Xt(ii,:,:));
        XCv=[XCv temp(:)];
    end
    [SE,tStat,pValue] = LinRegressionStat(beta,RS(:),k,XCv,0);

    betaCIC=beta;
    pValue_beta_Conflict(mm,:)=pValue(5:12);
    betaC(mm,:)=[beta(5:12)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Add Diesel
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    indx=ModelConflictDiesel{mm};
    % Test Dependence between Diesel and conflict
    [betap1,fval1,residual]=lsqnonlin(@(x)FitID(x,indx,betaM,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),log10(betaM(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
    
    [betap2,fval2,residual]=lsqnonlin(@(x)FitID(x,indx,betaCIC,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),log10(betaCIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
    if(fval1<fval2)
       betap=betap1;
    else
       betap=betap2; 
    end
    betap=10.^betap;
    beta=zeros(1,28);
    beta(indx)=betap;
    
    % Test for significance
    RS=WI(GNZI,(maxtau+1):end)-Yt;
    temp=squeeze(Xt(1,:,:));
    XCv=[temp(:)];
    for ii=2:32
        temp=squeeze(Xt(ii,:,:));
        XCv=[XCv temp(:)];
    end
    [SE,tStat,pValue] = LinRegressionStat(beta,RS(:),k,XCv,0);

    betaMIC=beta;
    pValue_beta_Mediation(mm,:)=pValue(5:12);
    
    betaC2(mm,:)=beta(5:12); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Relationship covariates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    W=Xt(17,:,:);
    W=W(:);
    if(ismember(5,ModelConflict{mm}))
        C=Xt(5,:,:);
        C=C(:);
        opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
        [bd,fval1]=fmincon(@(x)(mean((W-(10.^x(1)+10.^x(2).*C)).^2)),log10([162 32]),[],[],[],[],[-1 -1],[3 3],[],opts);

        RS=W-(10.^bd(1)+10.^bd(2).*C);
        temp=ones(size(C));
        XCv=[temp(:) C];    
        [SE,tStat,pValue_Corr(mm,1:2)] = LinRegressionStat(10.^bd,RS(:),2,XCv,0);
        betaCorr(mm,1:2)=10.^bd;
    end
    if(ismember(9,ModelConflict{mm}))
        C=Xt(9,:,:);
        C=C(:);
        opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
        [bd,fval1]=fmincon(@(x)(mean((W-(10.^x(1)+10.^x(2).*C)).^2)),log10([162 32]),[],[],[],[],[-1 -1],[3 3],[],opts);

        RS=W-(10.^bd(1)+10.^bd(2).*C);
        temp=ones(size(C));
        XCv=[temp(:) C];    
        [SE,tStat,pValue_Corr(mm,3:4)] = LinRegressionStat(10.^bd,RS(:),2,XCv,0);
        betaCorr(mm,3:4)=10.^bd;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Bootstrap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    betaCv=zeros(1000,8);
    betaC2v=zeros(1000,8);
    parfor ii=1:1000
        
        % conflict only
        GNZIt=GNZI(randi(length(GNZI),size(GNZI)));
        indx=ModelConflict{mm};
        % Test Dependence between Diesel and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaCIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),Temptv(GNZIt,:),r0,temp_0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaCIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,28);
        beta(indx)=betap;
    
        betaCv(ii,:)=beta(5:12); 
        
        
        indx=ModelConflictDiesel{mm};
        % Test Dependence between Diesel and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaMIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),Temptv(GNZIt,:),r0,temp_0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaMIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,32);
        beta(indx)=betap;
    
        betaC2v(ii,:)=beta(5:12); 
    end
    Boot_betaC{mm}=betaCv;
    Boot_betaC2{mm}=betaC2v;
    save('Mediation_Analysis_Alt_Wheat.mat');
end