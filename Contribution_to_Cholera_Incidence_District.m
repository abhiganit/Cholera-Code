function [CCRT,maxtau,GNZI,RC,MI] = Contribution_to_Cholera_Incidence_District()

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData;
% Loads the combinations of the different models
load('Combo.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load population size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('PopulationSize_DistrictYemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end); % Incidence data for the districts starts at may 1 2017

AICs=zeros(32,1);
Model_Inc=cell(32,1);
% Direct effect
CCR=cell(32,8);
MI_v=cell(32,1);
for kk=1:32
    load(['Fit-Vaccination-IncidenceperCapita' C(INC{kk}).N '-CalibratedDAR.mat'],'par','RSSv','XU');
    Model_Inc{kk}=[C(INC{kk}).N];
    
    [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,DAR,w,sigma_w]=RetParameterGA(par,XU,maxtau);
    
    n=n-TruncV;
    n(n<1)=1;
    
    [Yt,X]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
    dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
    dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
    EOVC=EffectOCV(dV1,KV,dV2,KV);

    MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

    for mm=1:8
        tempmat=zeros(size(squeeze(X(1,:,:))));
        for ii=(maxtau*(mm-1)+1):(mm.*maxtau)
            tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000;
        end
        CCR{kk,mm}=tempmat;
    end

    % Indirect effect

    if(ismember(4,INC{kk}) && sum(ismember([2 3],INC{kk}))>0) 
    %% Conflict indirect effect
        load(['Mediation_Analysis' C(INC{kk}).N '.mat'],'betaCorr');
        load(['Mediation_Analysis' C(INC{kk}).N '_District_Level.mat'],'C_MA');
        mmt=4;
        tempmat1=zeros(size(squeeze(X(1,:,:))));
        tempmat2=zeros(size(squeeze(X(1,:,:))));
        for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
            Diesel_Estimate=(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000;

            Y_Med_A_Tot=betaCorr(1).*ones(size(C_MA,2),size(C_MA,3));
            for jj=2:length(betaCorr)
                Y_Med_A_Tot=Y_Med_A_Tot+betaCorr(jj).*squeeze(C_MA(jj-1,:,:));
            end

            for zz=1:size(C_MA,1)                
                if(ismember(2,INC{kk}))
                    Conflict_1_Proportion=(betaCorr(zz+1).*squeeze(C_MA(zz,:,:)))./Y_Med_A_Tot;
                    tempmat1=tempmat1+Diesel_Estimate.*Conflict_1_Proportion;
                end
                if(ismember(3,INC{kk}))
                    Conflict_2_Proportion=(betaCorr(zz+1).*squeeze(C_MA(zz,:,:)))./Y_Med_A_Tot;
                    tempmat2=tempmat2+Diesel_Estimate.*Conflict_2_Proportion;
                end
            end
        end
        CCR{kk,2}=CCR{kk,2}+tempmat1;
        CCR{kk,3}=CCR{kk,3}+tempmat2;
        CCR{kk,4}=CCR{kk,4}-tempmat1-tempmat2;
    end
    
    if(ismember(5,INC{kk}) && sum(ismember([2 3],INC{kk}))>0) 
    %% Conflict indirect effect
        load(['Mediation_Analysis' C(INC{kk}).N '.mat'],'betaCorr');
        load(['Mediation_Analysis' C(INC{kk}).N '_District_Level.mat'],'C_MA');
        mmt=5;
        tempmat1=zeros(size(squeeze(X(1,:,:))));
        tempmat2=zeros(size(squeeze(X(1,:,:))));
        for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
            Wheat_Estimate=(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000;

            Y_Med_A_Tot=betaCorr(1).*ones(size(C_MA,2),size(C_MA,3));
            for jj=2:length(betaCorr)
                Y_Med_A_Tot=Y_Med_A_Tot+betaCorr(jj).*squeeze(C_MA(jj-1,:,:));
            end

            for zz=1:size(C_MA,1)                
                if(ismember(2,INC{kk}))
                    Conflict_1_Proportion=(betaCorr(zz+1).*squeeze(C_MA(zz,:,:)))./Y_Med_A_Tot;
                    tempmat1=tempmat1+Wheat_Estimate.*Conflict_1_Proportion;
                end
                if(ismember(3,INC{kk}))
                    Conflict_2_Proportion=(betaCorr(zz+1).*squeeze(C_MA(zz,:,:)))./Y_Med_A_Tot;
                    tempmat2=tempmat2+Wheat_Estimate.*Conflict_2_Proportion;
                end
            end
        end
        CCR{kk,2}=CCR{kk,2}+tempmat1;
        CCR{kk,3}=CCR{kk,3}+tempmat2;
        CCR{kk,5}=CCR{kk,5}-tempmat1-tempmat2;
    end
    
    
    MI_v{kk}=MI;
    AICs(kk)=aicbic(-RSSv,k);
end
AICs=AICs(1:32);
d_AIC=AICs-min(AICs);
w_AIC=exp(-d_AIC./2)./sum(exp(-d_AIC./2));

CCRT=cell(8,1);
MI=zeros(size(MI_v{1}));
for jj=1:8
    CCRT{jj}=CCR{1,jj}.*w_AIC(1);
    for mm=2:length(AICs)
        CCRT{jj}=CCRT{jj}+CCR{mm,jj}.*w_AIC(mm);
    end
end

for mm=1:32
   MI=MI+MI_v{mm}.*w_AIC(mm);
end

end

