clear;
clc;

FC=hex2rgb('#28595E');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenCountryData;

load('Fit-Vaccination-IncidenceperCapita-Targeted-Conflict-Diesel-Rain.mat');
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,DAR,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModelYemen(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_Gov_Incidence.mat')
IData=IData';
IData=sum(IData,1);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
AP=sum(AP,1);
PopS=[ repmat(AP(1),1,NW2016) repmat(AP(2),1,52)  repmat(AP(3),1,52)  repmat(AP(4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.065,0.62,0.93,0.35]);
NW=length(WI(1,:));
bar([(1+maxtau):NW],(MI),'Facecolor',FC,'LineStyle','none','Facealpha',0.35); hold on
scatter([1:NW],(IData),40,'k','filled'); 
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
IW=[1 21; 22 74; 75 121; 122 149];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',18);
ylim([0 55005]);
text(hy.Extent(1),56020,'A','Fontsize',32,'FontWeight','bold');


load('Fit-Vaccination-IncidenceperCapita-Targeted-Conflict-Diesel-Rain.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,DAR,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_Gov_Incidence.mat')
IData=IData';
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI2=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);


subplot('Position',[0.065,0.05,0.93,0.35]);
plot([(1+maxtau):NW],(MI)./sum(MI2)-1,'k','LineWidth',2); hold on
box off;
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTick',[-1:0.2:1],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
hy=ylabel([{'Error relative to'},{'aggregation of Governorates'}],'Fontsize',18);
xlabel('Week reported','Fontsize',18);
IW=[1 21; 22 74; 75 121; 122 149];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(-1,1,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),5020./5000,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),5020./5000,WN(4).N,'Fontsize',18);
ylim([-1 1]);
text(hy.Extent(1),5020./5000,'B','Fontsize',32,'FontWeight','bold');