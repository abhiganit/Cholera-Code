clear;
clc;

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:32];

AIC=zeros(length(INN),1);
Mean_Log_Likelihood_Validation=zeros(22,length(INN));
R2=zeros(22,length(INN)+1);
Targeted=zeros(length(INN),1);
Conflict=zeros(length(INN),1);
Shellings=zeros(length(INN),1);
Diesel=zeros(length(INN),1);
Wheat=zeros(length(INN),1);
k=zeros(length(INN),1);
MN=zeros(length(INN),1);

MN(end)=NaN;


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
District=cell(length(SD)+2,1);
Governorate=cell(length(SD)+2,1);
District(2:end-1)={SD.ADM2_EN};
Governorate(2:end-1)={SD.ADM1_EN};

District{end}='Hodeidah City';
Governorate{end}=SD(1).ADM1_EN;

Governorate{1}='AIC weight';
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
GNZI=GNZI(1:22);
NW=length(WI(1,:));
load('Combo.mat');


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

YtA=cell(length(INN),1);

for ii=1:length(INN)
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']))
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']);
        [lb,ub,IntC,pars] = BoundsFitting(XU,par,maxtau);
        
                
        [k(ii),beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,DAR,w,~]=RetParameterGA(par,XU,4);
        AIC(ii)=aicbic(-RSSv,k(ii));
        Log_Likelihood_Validation=-OFuncDistrict(pars,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),Temptv(GNZI,1:NW),PopS(GNZI,1:NW),CI(GNZI,1:NW));
        
        Mean_Log_Likelihood_Validation(:,ii)=Log_Likelihood_Validation./NW;
                
        n=n-TruncV;
        n(n<1)=1;
        [Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
        YtA{ii}=Yt;
        for mm=1:length(GNZI)
            R2(mm,ii)=corr(Yt(mm,:)',WI(GNZI(mm),(maxtau+1):end)').^2;
        end
    end
end


AIC=AIC-min(AIC);
wAIC=exp(-AIC./2)./sum(exp(-AIC./2));
wAIC=wAIC';

y_All=zeros(size(YtA{1}));
for jj=1:length(wAIC)
    y_All=y_All+wAIC(jj).*YtA{jj};
end
for mm=1:length(GNZI)
    R2(mm,end)=corr(y_All(mm,:)',WI(GNZI(mm),(maxtau+1):end)').^2;
end
        
Mean_Log_Likelihood_Validation=[wAIC;Mean_Log_Likelihood_Validation];
R2=[[wAIC 1];R2];
T=table(District,Governorate,Mean_Log_Likelihood_Validation);
T2=table(District,Governorate,R2);
writetable(T,'Model_Spatial_Validation_Log_likelihood.csv','Delimiter',',');
writetable(T2,'Model_Spatial_Validation_R2.csv','Delimiter',',');