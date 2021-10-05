close all;
clear;
FC=hex2rgb('#a6bddb');
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
% [WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDataVal;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

%% District

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

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
load('Yemen_District_Incidence.mat');
IData=IData';
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;
[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

SD=SD([29 31 71 fS' fA']);
for ii=1:(length(SD)) % Do not want to examine the last three as these are the larger areas and not the districts
    startDateofSim = datenum('5-01-2016');% Start date
    figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

    subplot('Position',[0.07,0.225,0.92,0.755]);
    dW=10;
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    NW=length(WI(1,:));
    bar([(1+maxtau):NW],MI(ii,:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
    scatter([1:NW],IData(ii,:),20,'k','filled'); 
    box off;
    xlim([0.5 length(WI(1,:))+0.5]);
    set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTick',[0:1000:8000],'YMinortick','on');
    ylim([0 max(ylim).*1.1])
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    xtickangle(90);
    hy=ylabel('Suspected cholera cases','Fontsize',18);
%     ylim([0 100]);

text(1.5,max(ylim),['\bf{\it{' SD(ii).ADM2_EN '}}'],'Fontsize',18);
    text(1.5,max(ylim)*(1-0.06),'Model estimate','Fontsize',18,'Color',FC(1,:));
%     text(1.5,max(ylim)*(1-0.12),'Model validation','Fontsize',18,'Color',FC(2,:));
    text(1.5,max(ylim)*(1-0.12),'Reported cases','Fontsize',18);
    

    xlabel('Week reported','Fontsize',18);
    
    print(gcf,['District_Validation_' SD(ii).ADM2_EN '.png'],'-dpng','-r600');
end
