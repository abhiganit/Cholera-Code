%% Runs the mediation analysis for diesel
clear;
close all;
clc;
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

load('Combo.mat');
CName={'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'};

% Models without Wheat but with Conflict or Shellings
ModelUW={};
ModelU={};
ModelConflictWheat={};
ModelConflict={};
count=0;
for ii=1:length(INC)
    test=INC{ii};
    if(ismember(5,test) & (ismember(2,test)| ismember(3,test)))
        count=count+1;
        ModelUW{count}=[test];
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
        ModelConflictWheat{count}=lt;
        ModelConflict{count}=ltc;
    end
end


for mm=1:length(ModelConflict)
    % Load the baseline model with Wheat and Conflict
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelUW{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [kD,betaM,~,~,~,~,~,KPw,~,~,~,~,~,~]=RetParameterGA(par,XU,maxtau);
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelU{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,sigma_w]=RetParameterGA(par,XU,maxtau);
    KP(2)=KPw(2);
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
    pValue_beta_Conflict=pValue(5:12);
    betaC=[beta(5:12)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Add Wheat
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    indx=ModelConflictWheat{mm};
    % Test Dependence between Wheat and conflict
    [betap1,fval1,residual]=lsqnonlin(@(x)FitID(x,indx,betaM,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),log10(betaM(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
    
    [betap2,fval2,residual]=lsqnonlin(@(x)FitID(x,indx,betaCIC,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w),log10(betaCIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
    if(fval1<fval2)
       betap=betap1;
    else
       betap=betap2; 
    end
    betap=10.^betap;
    beta=zeros(1,32);
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
    pValue_beta_Mediation=pValue(5:12);
    
    betaC2=beta(5:12); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Relationship covariates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    W=squeeze(Xt(17,:,:));
    
    
    XCF=[2 3];
    XCF=XCF(ismember([2 3],ModelU{mm}));
    
    C_MA=Xt(4.*(XCF-1)+1,:,:);
    
    opts=optimoptions('ga','MaxGenerations',10^4,'FunctionTolerance',10^(-9),'UseParallel',true); 
    
    [bd,fval1]=ga(@(x)Mediation_Fitting(x,W,C_MA),length(XCF)+1,[],[],[],[],-6.*ones(length(XCF),1),3.*ones(length(XCF),1),[],opts);
    
    beta_j=10.^bd;
    Y=beta_j(1).*ones(size(C_MA,2),size(C_MA,3));
    for jj=2:length(beta_j)
        Y=Y+beta_j(jj).*squeeze(C_MA(jj-1,:,:));
    end
    Y=Y(:);
    RS=(W(:)-Y(:));
    temp=ones(size(W(:)));
    
    XCv=[temp(:)];    
    for jj=1:length(XCF)
       temp_c=squeeze(C_MA(jj,:,:)); 
       XCv=[XCv temp_c(:)];    
    end
    [SE,tStat,pValue_Corr] = LinRegressionStat(10.^bd,RS(:),2,XCv,0);
    betaCorr=beta_j;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Bootstrap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    betaCv=zeros(1000,8);
    betaC2v=zeros(1000,8);
    parfor ii=1:1000
        
        % conflict only
        GNZIt=GNZI(randi(length(GNZI),size(GNZI)));
        indx=ModelConflict{mm};
        % Test Dependence between Wheat and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaCIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),Temptv(GNZIt,:),r0,temp_0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaCIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,32);
        beta(indx)=betap;
    
        betaCv(ii,:)=beta(5:12); 
        
        
        indx=ModelConflictWheat{mm};
        % Test Dependence between Wheat and conflict
        [betap,fval,residual]=lsqnonlin(@(x)FitID(x,indx,betaMIC,tA(GNZIt,:),DB,DA,Ctv(GNZIt,:),K,n,tau,maxtau,WPIN(GNZIt,:),FPIN(GNZIt,:),Mt(GNZIt,:),Wheatt(GNZIt,:),Dieselt(GNZIt,:),KP,V1(GNZIt,:),V2(GNZIt,:),KV,dV,Rtv(GNZIt,:),Temptv(GNZIt,:),r0,temp_0,WI(GNZIt,:),PopS(GNZIt,:),CI(GNZIt,:),DAR,w),log10(betaMIC(indx)),-32.*ones(size(indx)),3.*ones(size(indx)));
        betap=10.^betap;
        beta=zeros(1,32);
        beta(indx)=betap;
    
        betaC2v(ii,:)=beta(5:12); 
    end
    Boot_betaC=betaCv;
    Boot_betaC2=betaC2v;
    save(['Mediation_Analysis' CName{ModelUW{mm}}  '.mat'],'betaC','pValue_beta_Conflict','pValue_beta_Mediation','betaC2','betaCorr','pValue_Corr','C_MA','Boot_betaC','Boot_betaC2');
end