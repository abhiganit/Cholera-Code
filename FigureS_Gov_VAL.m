close all;
clear;
FC=[hex2rgb('#F5BE41');
    hex2rgb('#a6bddb');];
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDataVal;
NW=166; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
NW1=153;
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,~,w,sigma_w]=RetParameterGA(par,XU,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),Temptv(GNZI,:),r0,temp_0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_Gov_Incidence_Val.mat')
IData=IData';
IData=IData(GNZI,:);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=166-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

subplot('Position',[0.073684210526316,0.225,0.916315789473684,0.745443349753695]);
NW=length(WI(1,:));
b=bar([(1+maxtau):NW],sum(MI),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on

b.FaceColor = 'flat';
for mm=((NW1-maxtau)+1:162)
    b.CData(mm,:) = FC(2,:);
end

scatter([1:NW],sum(IData),40,'k','filled'); 
box off;
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTick',[0:5000:55000],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
hy=ylabel('Suspected cholera cases','Fontsize',18);
xlabel('Week reported','Fontsize',18);
IW=[1 21; 22 74; 75 121; 122 162];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',18);
ylim([0 55005]);
text(1.15,54005,'Model estimate','Fontsize',18,'Color',FC(1,:));
text(1.15,54005*(1-0.06),'Model validation','Fontsize',18,'Color',FC(2,:));
text(1.15,54005*(1-0.12),'Reported cases','Fontsize',18);
print(gcf,['Country_Validation_Yemen.png'],'-dpng','-r600');

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

S=S(GNZI);
for ii=1:length(GNZI)
    figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

    subplot('Position',[0.07,0.225,0.92,0.755]);
    dW=10;
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    NW=length(WI(1,:));
    b=bar([(1+maxtau):NW],MI(ii,:),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on

    b.FaceColor = 'flat';
    for mm=((NW1-maxtau)+1:162)
        b.CData(mm,:) = FC(2,:);
    end
    scatter([1:NW],IData(ii,:),20,'k','filled'); 
    box off;
    xlim([0.5 length(WI(1,:))+0.5]);
    
    set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','off','YMinortick','on');
    xlabel('Week reported','Fontsize',18);
    hy=ylabel('Suspected cholera cases','Fontsize',18);
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    xtickangle(90);
%     ylim([0 uub(ll)]);
    text(1.5,max(ylim),['\bf{\it{' S(ii).ADM1_EN '}}'],'Fontsize',18);
    text(1.5,max(ylim)*(1-0.06),'Model estimate','Fontsize',18,'Color',FC(1,:));
    text(1.5,max(ylim)*(1-0.12),'Model validation','Fontsize',18,'Color',FC(2,:));
    text(1.5,max(ylim)*(1-0.18),'Reported cases','Fontsize',18);
    print(gcf,['Gov_Validation_' S(ii).ADM1_EN '.png'],'-dpng','-r600');
end
