%% Runs the mediation analysis for diesel
clear;
close all;
clc;
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData;

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
    if(ismember(4,test) & (ismember(2,test)| ismember(3,test)))
        count=count+1;
        ModelUD{count}=[test];
        lt=[];
        ltc=[];
        tt=[];
        for jj=1:length(test)
            lt=[lt 4.*(test(jj)-1)+[1:4]];
            if(test(jj)~=4)
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


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;


for mm=1:length(ModelConflict)
    % Load the baseline model with Diesel and Conflict
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelUD{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [kD,betaM,~,~,~,~,~,KPd,~,~,~,~,~,~]=RetParameterGA(par,XU,maxtau);
    load(['Fit-Vaccination-IncidenceperCapita' CName{ModelU{mm}} '-CalibratedDAR.mat'],'par','XU','DAR')    
    [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,sigma_w]=RetParameterGA(par,XU,maxtau);
    KP(1)=KPd(1);
    
    n=n-TruncV;
    n(n<1)=1;
    [Yt,Xt]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
        
    XCF=[2 3];
    test=ismember([2 3],ModelU{mm});
    XCF=XCF(test);
    C_MA=Xt(4.*(XCF-1)+1,:,:);
    save(['Mediation_Analysis' CName{ModelUD{mm}}  '_District_Level.mat'],'C_MA'); %,'Boot_betaC','Boot_betaC2'
end