%% AIC MSE CVE for the various models
clear;
clc;
[WIF,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZIF,GV,maxtau,PopS,CI] = LoadYemenData;
NW=length(WIF(1,:));

[WIVal,CtvVal,tAVal,RtvVal,MtVal,PVal,RCVal,HVal,WPINVal,FPINVal,DieseltVal,WheattVal,V1Val,V2Val,GNZIVal,GVVal,maxtauVal,PopSVal,CIVal] = LoadYemenDataVal;
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
close all;
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:64];
nd=WIF(GNZIF,(maxtau+1):end);
nd=length(nd(:));
CVES=10^6.*ones(length(INN)+1,1);
CVET=10^6.*ones(length(INN)+1,1);
BIC=10^6.*ones(length(INN)+1,1);
AIC=10^6.*ones(length(INN)+1,1);
MSE=10^6.*ones(length(INN)+1,1);
MeanDataFit=10^6.*ones(length(INN)+1,1);
MeanDataCVS=10^6.*ones(length(INN)+1,1);
MeanDataCVT=10^6.*ones(length(INN)+1,1);
Targeted=zeros(length(INN)+1,1);
Conflict=zeros(length(INN)+1,1);
Shellings=zeros(length(INN)+1,1);
Diesel=zeros(length(INN)+1,1);
Wheat=zeros(length(INN)+1,1);
Rain=zeros(length(INN)+1,1);
k=zeros(length(INN)+1,1);
MN=zeros(length(INN)+1,1);

MN(end)=NaN;
MSE(end)=NaN;
CVES(end)=NaN;
CVET(end)=NaN;


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
NW=length(WI(1,:));
load('Combo.mat');
for ii=1:length(INN)
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']))
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']);
        [lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par,CF,maxtau);
        MSE(ii)=RSSv;
        
        tempD=WIF(GNZIF,(maxtau+1):end);
        MeanDataFit(ii)=mean(tempD(:));
        
        tempD=WI(GNZI,(maxtau+1):end);        
        MeanDataCVS(ii)=mean(tempD(:));
        
        tempD=WIVal(GNZIVal,(154):end);        
        MeanDataCVT(ii)=mean(tempD(:));
        
        [k(ii)]=RetParameterPS(par,XU,CF,4);
        AIC(ii)=AICScore(k(ii),nd,RSSv);
%         BIC(ii)=BICScore(k(ii),nd,RSSv);
        CVES(ii)=OFuncDistrict(pars,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW));
        CVET(ii)=OFuncProGA_TempVal(pars,CF,WIVal(GNZIVal,1:end),tAVal(GNZIVal,1:end),CtvVal(GNZIVal,1:end),XU,maxtau,WPINVal(GNZIVal,1:end),FPINVal(GNZIVal,1:end),MtVal(GNZIVal,1:end),WheattVal(GNZIVal,1:end),DieseltVal(GNZIVal,1:end),V1Val(GNZIVal,1:end),V2Val(GNZIVal,1:end),RtvVal(GNZIVal,1:end),RF,PopSVal(GNZIVal,1:end),CIVal(GNZIVal,1:end));
        if(ismember(1,INC{INN(ii)}))
            Targeted(ii)=1;
        else
            Targeted(ii)=0;
        end
        if(ismember(2,INC{INN(ii)}))
            Conflict(ii)=1;
        else
            Conflict(ii)=0;
        end
        if(ismember(3,INC{INN(ii)}))
            Shellings(ii)=1;
        else
            Shellings(ii)=0;
        end
        if(ismember(4,INC{INN(ii)}))
            Diesel(ii)=1;
        else
            Diesel(ii)=0;
        end
        if(ismember(5,INC{INN(ii)}))
            Wheat(ii)=1;
        else
            Wheat(ii)=0;
        end
        if(ismember(6,INC{INN(ii)}))
            Rain(ii)=1;
        else
            Rain(ii)=0;
        end
        MN(ii)=ii;
    end
end
 AIC(1:end-1)=AIC(1:end-1)-min(AIC(1:end-1));
 AIC(end)=NaN;
 wAIC=exp(-AIC(1:end-1)./2)./sum(exp(-AIC(1:end-1)./2));
 wAIC=[wAIC;1];
%  BIC=BIC-min(BIC);
%  wBIC=exp(-BIC./2)./sum(exp(-BIC./2));
 
 NormMSE=MSE./MeanDataFit; 
 NormCVS=CVES./MeanDataCVS;
 NormCVT=CVET./MeanDataCVT;
Targeted(end)=sum(Targeted(1:end-1).*wAIC(1:end-1));
Conflict(end)=sum(Conflict(1:end-1).*wAIC(1:end-1));
Shellings(end)=sum(Shellings(1:end-1).*wAIC(1:end-1));
Diesel(end)=sum(Diesel(1:end-1).*wAIC(1:end-1));
Wheat(end)=sum(Wheat(1:end-1).*wAIC(1:end-1));
Rain(end)=sum(Rain(1:end-1).*wAIC(1:end-1));

 T=table(MN,Targeted,Conflict,Shellings,Diesel,Wheat,Rain,MSE,CVES,CVET,NormMSE,NormCVS,NormCVT,k,AIC,wAIC);
 
 writetable(T,'ModelFit.csv','Delimiter',',');
 
 BestM=find(AIC==0);
 
 clearvars -except BestM INC INN C
 
 load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(BestM)}).N '-CalibratedDAR.mat']);
 
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
 % Governoerate 
 S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
S=S(GNZI);
% Record the names for the IDP calculation
Sm=cell(length(S),1); % allocate space
for ii=1:length(Sm)
Sm{ii}=S(ii).ADM1_EN; % record name
end


SD = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen

fS=zeros(length(SD),1);
for ii=1:length(fS)
  fS(ii)=strcmp({'Amanat Al Asimah'},SD(ii).ADM1_EN); 
end
fS=find(fS==1);

fA=zeros(length(SD),1);
for ii=1:length(fA)
  fA(ii)=strcmp({'Aden'},SD(ii).ADM1_EN); 
end

fA=find(fA==1);


SD=SD([29 31 71 fS' fA']);

MSE=zeros(length(Sm)+length(SD)+1,1);
Mean_Data=zeros(length(Sm)+length(SD)+1,1);
R2=zeros(length(Sm)+length(SD)+1,1);


NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,DAR,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

MSE(1:length(Sm))=mean((Yt-WI(GNZI,(maxtau+1):end)).^2,2);
Mean_Data(1:length(Sm))=mean(WI(GNZI,(maxtau+1):end),2);
for ii=1:length(Sm)
    R2(ii)=corr(Yt(ii,:)',WI(GNZI(ii),(maxtau+1):end)').^2;
end


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,DAR,w]=RetParameterPS(par,XU,CF,maxtau);


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;

[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

MSE((1+length(Sm)):end)=mean((Yt(1:end-2,:)-WI(GNZI(1:end-2,:),(maxtau+1):end)).^2,2);
Mean_Data((1+length(Sm)):end)=mean(WI(GNZI(1:end-2),:),2);

for ii=(1+length(Sm)):length(MSE)
    R2(ii)=corr(Yt(ii-length(Sm),:)',WI(GNZI(ii-length(Sm)),(maxtau+1):end)').^2;
end

Norm_MSE=MSE./Mean_Data;
Location={Sm{:},SD.ADM2_EN,'Hodeidah City'}';


 T=table(Location,MSE,Mean_Data,Norm_MSE,R2);
 
 writetable(T,'AICModelStat.csv','Delimiter',',');
        