function [Y_Avg,NW,NW1,IData_VAL,maxtau,GNZI] = Weighted_Model_Incidence(Level_Fit)

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:32];

MIv=cell(1,32);
AIC_s=zeros(1,32);
load('Combo.mat');

if(strcmp(Level_Fit,'Governorate'))
    load('PopulationSize_Yemen.mat');
    NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
    NW2019=166-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
    % External effect due to IDP
    PopS_Val=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

    [WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDataVal;
    NW=length(WI(1,:));
    NW1=153;
    load('Yemen_Gov_Incidence_Val.mat')
    IData_VAL=IData';
    IData_VAL=IData_VAL(GNZI,:);
    for mm=1:32
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(mm)}).N '-CalibratedDAR.mat']);
        
        % Evaluate the number of paramters that are being used in the estimation 
        [k_par,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,~]=RetParameterGA(par,XU,maxtau);
        AIC_s(mm)=aicbic(-RSSv,k_par);

        [Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
        MIv{mm}=(Yt./(10000)).*PopS_Val(GNZI,maxtau+1:end);
    end
    
    d_AIC=AIC_s-min(AIC_s);
    w=exp(-d_AIC./2)./sum(exp(-d_AIC./2));
    Y_Avg=zeros(size(MIv{1}));
    for mm=1:32
        Y_Avg=Y_Avg+w(mm).*MIv{mm};
    end
elseif(strcmp(Level_Fit,'District'))
    
    [WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
    load('Yemen_District_Incidence.mat');
    IData_VAL=IData';
    NW=length(WI(1,:));
    NW1=NW;
    startDateofSim = datenum('10-03-2016');% Start date
    endDateofSim = datenum('5-01-2017');% End date
    TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

    
    for mm=1:32
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(mm)}).N '-CalibratedDAR.mat']);
        
        % Evaluate the number of paramters that are being used in the estimation 
        [k_par,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,~]=RetParameterGA(par,XU,maxtau);
        AIC_s(mm)=aicbic(-RSSv,k_par);
        
        
        n=n-TruncV;
        n(n<1)=1;
        [Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
        MIv{mm}=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
    end
    
    d_AIC=AIC_s-min(AIC_s);
    w=exp(-d_AIC./2)./sum(exp(-d_AIC./2));
    Y_Avg=zeros(size(MIv{1}));
    for mm=1:32
        Y_Avg=Y_Avg+w(mm).*MIv{mm};
    end
    
end

end

