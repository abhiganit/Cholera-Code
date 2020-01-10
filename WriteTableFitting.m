%% AIC MSE CVE for the various models
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=length(WI(1,:));
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
close all;
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:64];
nd=WI(GNZI,(maxtau+1):end);
nd=length(nd(:));
CVE=10^6.*ones(length(INN),1);
AIC=10^6.*ones(length(INN),1);
MSE=10^6.*ones(length(INN),1);
Targeted=zeros(length(INN),1);
Conflict=zeros(length(INN),1);
Shellings=zeros(length(INN),1);
Diesel=zeros(length(INN),1);
Wheat=zeros(length(INN),1);
Rain=zeros(length(INN),1);
k=zeros(length(INN),1);

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
NW=length(WI(1,:));
load('Combo.mat');
for ii=1:length(INN)
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']))
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']);
        [lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par,CF,maxtau);
        MSE(ii)=RSSv;
        [k(ii)]=RetParameterPS(par,XU,CF,4);
        AIC(ii)=AICScore(k(ii),nd,RSSv);
        CVE(ii)=OFuncDistrict(pars,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW));
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
    end
end
 AIC=AIC-min(AIC);

 T=table(Targeted,Conflict,Shellings,Diesel,Wheat,Rain,MSE,CVE,k,AIC);
 
 writetable(T,'ModelFit.csv','Delimiter',',');
 
 BestM=find(AIC==0);
 
 clearvars -except BestM INC INN C
 
 load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(BestM)}).N '.mat']);
 
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
R2=zeros(length(Sm)+length(SD)+1,1);


NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,DAR,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

MSE(1:length(Sm))=mean((Yt-WI(GNZI,(maxtau+1):end)).^2,2);
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
for ii=(1+length(Sm)):length(MSE)
    R2(ii)=corr(Yt(ii-length(Sm),:)',WI(GNZI(ii-length(Sm)),(maxtau+1):end)').^2;
end


Location={Sm{:},SD.ADM2_EN,'Hodeidah City'}';


 T=table(Location,MSE,R2);
 
 writetable(T,'AICModelStat.csv','Delimiter',',');
        