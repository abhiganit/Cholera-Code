close all;
clear;
FC=hex2rgb('#28595E');
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);
% 
% 
% [Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
% load('Yemen_Gov_Incidence.mat')
% % IData=IData';
% % load('PopulationSize_Yemen.mat');
% % NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
% % NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% % % External effect due to IDP
% % PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
% 
% MI=(Yt);
% 
% S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
% S=S(GNZI);
% IData=IData(GNZI,:);
% WI=WI(GNZI,:);
% for ii=1:11
%     startDateofSim = datenum('10-03-2016');% Start date
%     figure('units','normalized','outerposition',[0 0 1 1]);
% 
%     subplot('Position',[0.045,0.62,0.46,0.35]);
%     dW=10;
%     XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
%     NW=length(WI(1,:));
%     bar([(1+maxtau):NW],MI(1+2.*(ii-1),:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
%     scatter([1:NW],WI(1+2.*(ii-1),:),20,'k','filled'); 
%     box off;
%     xlim([0.5 length(WI(1,:))+0.5]);
%     set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTick',[0:1000:8000],'YMinortick','on');
%     % Sets the y-axis to not have 10^n
%     ax=gca; % finds the current axis
%     ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
%     xtickangle(90);
%     hy=ylabel('Incidence per captia','Fontsize',18);
%     ylim([0 max(WI(:))]);
%     text(1.5,max(WI(:)),S(1+2.*(ii-1)).ADM1_EN,'Fontsize',18);
% 
%     xlabel('Week reported','Fontsize',18);
% 
%     if(ii<11)
%         
%         subplot('Position',[0.535,0.62,0.46,0.35]);
%         NW=length(WI(1,:));
%         bar([(1+maxtau):NW],MI(2+2.*(ii-1),:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
%         scatter([1:NW],WI(2+2.*(ii-1),:),20,'k','filled'); 
%         box off;
%         xlim([0.5 length(WI(1,:))+0.5]);
%         set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTickLabel','','YTick',[0:1000:8000],'YMinortick','on');
%         % Sets the y-axis to not have 10^n
%         ax=gca; % finds the current axis
%         ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
%         xtickangle(90);
%         ylim([0 max(WI(:))]);
%         text(1.5,max(WI(:)),S(2+2.*(ii-1)).ADM1_EN,'Fontsize',18);
% 
%         xlabel('Week reported','Fontsize',18);
%     end
% end

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

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;
[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
MI=(Yt);
SD=SD([29 31 71 fS' fA']);
for ii=1:11
    startDateofSim = datenum('5-01-2016');% Start date
    figure('units','normalized','outerposition',[0 0 1 1]);

    subplot('Position',[0.045,0.62,0.46,0.35]);
    dW=10;
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    NW=length(WI(1,:));
    bar([(1+maxtau):NW],MI(1+2.*(ii-1),:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
    scatter([1:NW],WI(1+2.*(ii-1),:),20,'k','filled'); 
    box off;
    xlim([0.5 length(WI(1,:))+0.5]);
    set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTick',[0:1000:8000],'YMinortick','on');
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    xtickangle(90);
    hy=ylabel('Incidence per captia','Fontsize',18);
    ylim([0 100]);
    text(1.5,100,SD(1+2.*(ii-1)).ADM2_EN,'Fontsize',18);

    xlabel('Week reported','Fontsize',18);

    if(ii<11)
        
        subplot('Position',[0.535,0.62,0.46,0.35]);
        NW=length(WI(1,:));
        bar([(1+maxtau):NW],MI(2+2.*(ii-1),:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
        scatter([1:NW],WI(2+2.*(ii-1),:),20,'k','filled'); 
        box off;
        xlim([0.5 length(WI(1,:))+0.5]);
        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTickLabel','','YTick',[0:1000:8000],'YMinortick','on');
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        xtickangle(90);
        ylim([0 100]);
        text(1.5,100,SD(2+2.*(ii-1)).ADM2_EN,'Fontsize',18);

        xlabel('Week reported','Fontsize',18);
    else
        subplot('Position',[0.535,0.62,0.46,0.35]);
        NW=length(WI(1,:));
        bar([(1+maxtau):NW],MI(2+2.*(ii-1),:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
        scatter([1:NW],WI(2+2.*(ii-1),:),20,'k','filled'); 
        box off;
        xlim([0.5 length(WI(1,:))+0.5]);
        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTickLabel','','YTick',[0:1000:8000],'YMinortick','on');
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        xtickangle(90);
        ylim([0 100]);
        text(1.5,100,'Hodeihdah City','Fontsize',18);

        xlabel('Week reported','Fontsize',18);
    end
end
