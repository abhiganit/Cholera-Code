close all;
clear;
[WI2,~,tA2,Rtv2,~,P2,RC2,H2,WPINm2,FPINm2,Dieselt2,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
RC=RC2;
WI2=WI2(GNZI,:);
tA2=tA2(GNZI,:);
Rtv2=Rtv2(GNZI,:);
P2=P2(GNZI,:);
RC2=RC2(GNZI);
H2=H2(GNZI);
WPINm2=WPINm2(GNZI,:);
FPINm2=FPINm2(GNZI,:);
Dieselt2=Dieselt2(GNZI,:);
Wheatt2=Wheatt2(GNZI,:);
% Need to increase the diesel price as we transformed it by subtracting the
% minimum
load('Diesel_Gov_Yemen.mat')
Dieselt2=Dieselt2+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt2=Wheatt2+min(Wheat(Wheat>0));

% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv2=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt2=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

Ctv2=Ctv2(GNZI,:);
Mt2=Mt2(GNZI,:);


S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


%% Incidence
load('Yemen_Gov_Incidence.mat')
IData=IData';
IData=IData(GNZI,:);
NW=153;
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(GNZI,:);
PopH=sum(PopS(RC2==1,:),1);
PopG=sum(PopS(RC2==0,:),1);
TestH=10000.*sum(IData(RC2==1,:),1)./PopH;
TestG=10000.*sum(IData(RC2==0,:),1)./PopG;

figure('units','normalized','outerposition',[0.1 0 0.55 1]);
subplot('Position',[0.065,0.863,0.43,0.13]); 
plot([1:153],TestH,'color',hex2rgb('#C60000'),'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#C60000'),'MarkerEdgeColor',hex2rgb('#C60000'), 'MarkerIndices', 1:2:153); hold on;
plot([1:153],TestG,'color',hex2rgb('#C60000'),'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#C60000'), 'MarkerIndices', 1:2:153);
startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
yh=ylabel({'Incidence per 10,000'},'Fontsize',10);
legend({'Houthi control','Government control'},'Fontsize',8)
legend boxoff;
text(-0.125,0.99,'A','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Rainfall
TestH=mean(Rtv2(RC2==1,:),1);
TestG=mean(Rtv2(RC2==0,:),1);

subplot('Position',[0.56,0.863,0.43,0.13]); 

plot([1:153],TestH,'color',[5,112,176]./255,'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',[5,112,176]./255,'MarkerEdgeColor',[5,112,176]./255, 'MarkerIndices', 1:2:153); hold on;
plot([1:153],TestG,'color',[5,112,176]./255,'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',[5,112,176]./255, 'MarkerIndices', 1:2:153);

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
box off;
xlabel('Week of report','Fontsize',10);
ylabel('Rainfall (mm)','Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8)
legend boxoff;
text(-0.125,0.99,'B','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Conflict


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Ctv2(RC2==1,:),1);
TestG=mean(Ctv2(RC2==0,:),1);

subplot('Position',[0.065,0.667,0.43,0.13]);  

plot([1:153],TestH,'color',hex2rgb('#DE7A22'),'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#DE7A22'),'MarkerEdgeColor',hex2rgb('#DE7A22'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#DE7A22'),'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#DE7A22'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel({'Number of','conflict events'},'Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'C','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Shelling

subplot('Position',[0.56,0.667,0.43,0.13]);  

TestH=mean(Mt2(RC2==1,:),1);
TestG=mean(Mt2(RC2==0,:),1);

plot([1:153],TestH,'color',hex2rgb('#4C3F54'),'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#4C3F54'),'MarkerEdgeColor',hex2rgb('#4C3F54'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#4C3F54'),'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#4C3F54'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel({'Number of','shelling attacks'},'Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'D','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Diesel
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Dieselt2(RC2==1,:),1);
TestG=mean(Dieselt2(RC2==0,:),1);

subplot('Position',[0.065,0.47,0.43,0.13]);  

plot([1:153],TestH,'color',[153,52,4]./255,'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',[153,52,4]./255,'MarkerEdgeColor',[153,52,4]./255, 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',[153,52,4]./255,'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',[153,52,4]./255, 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel('Diesel price','Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'E','Fontsize',18,'FontWeight','bold','Units','Normalized');
%% Wheat

subplot('Position',[0.56,0.47,0.43,0.13]);  

TestH=mean(Wheatt2(RC2==1,:),1);
TestG=mean(Wheatt2(RC2==0,:),1);

plot([1:153],TestH,'color',hex2rgb('#FAAF08'),'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#FAAF08'),'MarkerEdgeColor',hex2rgb('#FAAF08'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#FAAF08'),'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#FAAF08'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel('Wheat price','Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'F','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% target attacks


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019

subplot('Position',[0.065,0.27,0.43,0.13]);  

TestH=mean(tA2(RC2==1,:),1);
TestG=mean(tA2(RC2==0,:),1);
b=bar([1:153],[TestH; TestG]','stacked');

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel({'Number of attacks','on water systems'},'Fontsize',10);

ColorM=[[221,28,119]./255; [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=[221,28,119]./255;
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':'; 
legend('Houthi control','Government control');
legend boxoff;

text(-0.125,0.99,'G','Fontsize',18,'FontWeight','bold','Units','Normalized');
%% WaSH
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP

subplot('Position',[0.56,0.27,0.43,0.13]); 
load('WASH_PIN_Yemen.mat')
WPIN=WPIN(GNZI,:);
AP=AP(GNZI,:);
TestH=sum(WPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
TestG=sum(WPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
b=bar([2016 2017 2018 2019],[TestH;TestG]');

xlim([2015.5 2019.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',8,'YMinortick','on');
xtickangle(45);
xlabel('Year','Fontsize',10);
box off;
ylabel({'Proportion in','need of WaSH'},'Fontsize',10);

ColorM=[[0 0.6 0.6]; [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=[0 0.6 0.6];
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':';
ylim([0 0.8]);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'H','Fontsize',18,'FontWeight','bold','Units','Normalized');

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
%% Food security

subplot('Position',[0.065,0.07,0.43,0.13]); 
load('Food_PIN_Yemen.mat')
FPIN=FPIN(GNZI,:);
AP=AP(GNZI,:);
TestH=sum(FPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
TestG=sum(FPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
b=bar([2016 2017 2018 2019],[TestH;TestG]');

xlim([2015.5 2019.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',8,'YMinortick','on');
xtickangle(45);
xlabel('Year','Fontsize',10);
box off;
ylabel({'Proportion in need','of food security'},'Fontsize',10);

ColorM=[hex2rgb('#2E4600'); [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=hex2rgb('#2E4600');
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':';


legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
text(-0.125,0.99,'I','Fontsize',18,'FontWeight','bold','Units','Normalized');
%% Vaccination
[V1,V2]= VaccinationTime(1,153);
V1=V1(GNZI,:); % Received at least one dose
TempH=cumsum(sum(V1(RC2==1,:),1))./PopH;
TempG=cumsum(sum(V1(RC2==0,:),1))./PopG;


subplot('Position',[0.56,0.07,0.43,0.13]); 

plot([1:153],100.*TempH,'color',hex2rgb('#80BD9E'),'LineWidth',2,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#80BD9E'),'MarkerEdgeColor',hex2rgb('#80BD9E'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],100.*TempG,'color',hex2rgb('#80BD9E'),'LineWidth',1.75,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#80BD9E'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',8,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',10);
box off;
ylabel('Vaccination uptake','Fontsize',10);

legend({'Houthi control','Government control'},'Fontsize',8,'Location','NorthWest')
legend boxoff;
ytickformat('percentage');
text(-0.125,0.99,'J','Fontsize',18,'FontWeight','bold','Units','Normalized');
print(gcf,['Figure1'],'-depsc','-r600');
print(gcf,['Figure1'],'-dpng','-r600');