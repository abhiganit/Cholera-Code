%% Runs the mediation analysis for diesel
clear;
close all;
clc;
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData;

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


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

for mm=1:length(ModelConflict)
    % Load the baseline model with Wheat and Conflict
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelUW{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [kD,betaM,~,~,~,~,~,KPw,~,~,~,~,~,~]=RetParameterGA(par,XU,maxtau);
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelU{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,sigma_w]=RetParameterGA(par,XU,maxtau);
    KP(2)=KPw(2);
    
    n=n-TruncV;
    n(n<1)=1;
    [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
   
    XCF=[2 3];
    XCF=XCF(ismember([2 3],ModelU{mm}));
    
    C_MA=Xt(4.*(XCF-1)+1,:,:);
    
    save(['Mediation_Analysis' CName{ModelUW{mm}}  '_District_Level.mat'],'C_MA'); %,'Boot_betaC','Boot_betaC2'
end