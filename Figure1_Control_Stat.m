close all;
clear;
[WI2,~,tA2,Rtv2,~,P2,RC2,H2,WPINm2,FPINm2,Dieselt2,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
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

[WI,~,tA,Rtv,~,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
% Need to increase the diesel price as we transformed it by subtracting the
% minimum
load('Diesel_Gov_Yemen.mat')
Dieselt=Dieselt+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt=Wheatt+min(Wheat(Wheat>0));

% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

figure('units','normalized','outerposition',[0 0 1 1]);
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

TestH=WI2(RC2==1,:);
TestG=WI2(RC==0,:);
[h,pInc]=ttest2(TestH,TestG);
subplot(4,2,1);plot([1:153],pInc,'color',hex2rgb('#C60000'),'LineWidth',2); ylim([0 1]); box off;
%% Rainfall
TestH=(Rtv2(RC2==1,:));
TestG=(Rtv2(RC2==0,:));
[h,pRF]=ttest2(TestH,TestG);
subplot(4,2,2);plot([1:153],pRF,'color',[5,112,176]./255,'LineWidth',2); ylim([0 1]); box off;
%% Conflict

TestH=(Ctv2(RC2==1,:));
TestG=(Ctv2(RC2==0,:));
[h,pC]=ttest2(TestH,TestG);
subplot(4,2,3);plot([1:153],pC,'color',hex2rgb('#DE7A22'),'LineWidth',2); ylim([0 1]); box off;
%% Shelling

TestH=(Mt2(RC2==1,:));
TestG=(Mt2(RC2==0,:));

[h,pS]=ttest2(TestH,TestG);
subplot(4,2,4);plot([1:153],pS,'color',hex2rgb('#4C3F54'),'LineWidth',2); ylim([0 1]); box off;
%% Diesel
TestH=(Dieselt2(RC2==1,:));
TestG=(Dieselt2(RC2==0,:));

[h,pD]=ttest2(TestH,TestG);
subplot(4,2,5);plot([1:153],pD,'color',[153,52,4]./255,'LineWidth',2); ylim([0 1]); box off;
%% Wheat

TestH=(Wheatt2(RC2==1,:));
TestG=(Wheatt2(RC2==0,:));

[h,pW]=ttest2(TestH,TestG);
subplot(4,2,6);plot([1:153],pW,'color',hex2rgb('#FAAF08'),'LineWidth',2); ylim([0 1]); box off;
%% target attacks

TestH=(tA2(RC2==1,:));
TestG=(tA2(RC2==0,:));

[h,pT]=ttest2(TestH,TestG);
subplot(4,2,7);plot([1:153],pT,'color',[221,28,119]./255,'LineWidth',2); ylim([0 1]); box off;
%% WaSH
load('PopulationSize_Yemen.mat');
load('WASH_PIN_Yemen.mat')
WPIN=WPIN(GNZI,:);
AP=AP(GNZI,:);
TestH=(WPIN(RC2==1,:))./(AP(RC2==1,:));
TestG=(WPIN(RC2==0,:))./(AP(RC2==0,:));

[h,pWSH]=ttest2(TestH,TestG);
%% Food security
load('Food_PIN_Yemen.mat')
FPIN=FPIN(GNZI,:);
TestH=(FPIN(RC2==1,:))./(AP(RC2==1,:));
TestG=(FPIN(RC2==0,:))./(AP(RC2==0,:));

[h,pF]=ttest2(TestH,TestG);