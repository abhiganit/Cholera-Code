%% Inlcudes the effwects of conflict and shellngs on the diesel prices
% Read table of past fitsclose all;
close all;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
S=S(GNZI);

[Yt,X]= LogisticModelR(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(WI(GNZI,:)./(10000)).*PopS(GNZI,:);


figure('units','normalized','outerposition',[0 0 1 1]);

for jj=1:7
    for ii=1:3
        subplot('Position',[0.04+(0.33*(ii-1)),0.87-(0.16*(jj-1)),0.29,0.12]);
        plot(Yt(ii+3*(jj-1),:),'k','LineWidth',2)
        hold on
        plot(ones(size(Yt(ii+3*(jj-1),:))),'r-.','LineWidth',1);
        xlabel('time');
        ylabel('R_E');
        box off;
        set(gca,'tickdir','out','linewidth',2,'Fontsize',16);
        ylim([0 2]);
        xlim([1 149]);
        text(2.5,1.9,S(ii+3*(jj-1)).ADM1_EN,'Fontsize',14);
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);

for jj=1:7
    for ii=1:3
        subplot('Position',[0.04+(0.33*(ii-1)),0.87-(0.16*(jj-1)),0.29,0.12]);
        plot(MI(ii+3*(jj-1),(maxtau+1):153),'b','LineWidth',2)
        xlabel('time');
        ylabel('Cases');
        box off;
        set(gca,'tickdir','out','linewidth',2,'Fontsize',16,'Ytick',[]);
        %ylim([0 2]);
        xlim([1 149]);
        text(2.5,max(ylim)*1.9/2,S(ii+3*(jj-1)).ADM1_EN,'Fontsize',14);
    end
end