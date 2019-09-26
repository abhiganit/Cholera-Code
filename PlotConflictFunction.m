close all;
clear;
clc;
load('Yemen_Gov_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

WI=IData'; % Transpose the data set such that the number of areas is the row
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection

AF=2;
XU=ones(1,11);
load([pwd '\Tables\BackSelectModel-PercentDataSet=100-XU=2047-CF=1-RIF=0-PF=1-tau=1  1  1  3  1  4  4  1  2  1.mat'],'par');


%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
 
figure('units','normalized','outerposition',[0 0 1 1]);
    
subplot('Position',[0.042,0.163120567375887,0.29,0.4]);

R=ImpactRainfall(linspace(0,3,10001),RF,rh);

plot(linspace(0,3,10001),R,'color',hex2rgb('9ecae1'),'LineWidth',2); hold on;
    R=ImpactRainfall(linspace(0,1,10001),RIF,rl);
plot(linspace(0,3,10001),R,'color',hex2rgb('2171b5'),'LineWidth',2); hold off;

yh=ylabel('Rainfall function','Fontsize',18);
xlabel('Average rainfall per week','Fontsize',18);

legend({'Rainfall','Rainfall and Incidence'},'location','northwest');
legend boxoff;
box off;
set(gca,'Fontsize',16,'tickdir','out','LineWidth',2,'XTick',[0:0.2:3],'Xminortick','on','YTick',[0:0.2:3],'Yminortick','on');
text(yh.Extent(1),0.98*max(ylim),'A','Fontsize',32,'FontWeight','bold');
   C=ImpactConflict([0:75],K,n,CF);

subplot('Position',[0.3696,0.163120567375887,0.29,0.4]);
plot([0:75],C,'color',hex2rgb('67000d'),'LineWidth',2)
yh=ylabel('Conflit function','Fontsize',18);
xlabel('Number of conflict events per week','Fontsize',18);
set(gca,'Fontsize',16,'tickdir','out','LineWidth',2,'XTick',[0:5:75],'Xminortick','on','YTick',[0:2:20],'Yminortick','on');
box off;
ylim([0 20])
xlim([0 20])
text(yh.Extent(1),0.98*max(ylim),'B','Fontsize',32,'FontWeight','bold');
A=ImpactAttack([zeros(1,21) 1 zeros(1,21)],0,DAE,0,0);
At=ImpactAttack([zeros(1,21) 1 zeros(1,21)],DB,DA,0,0);
subplot('Position',[0.7,0.163120567375887,0.29,0.4]);
h=bar([-21:21],[A; At]','LineStyle','none'); 
xlim([-3.5 6.5]);
ylim([0 1.01])
h(2).FaceColor = 'flat';
    h(2).CData = hex2rgb('ef3b2c');
    h(1).FaceColor = 'flat';
    h(1).CData = hex2rgb('fc9272');
yh=ylabel('Attack function','Fontsize',18);
xlabel('Weeks from an attack','Fontsize',18);
legend({'Attack','Attack and Incidence'});
legend boxoff;
box off;
set(gca,'Fontsize',16,'tickdir','out','LineWidth',2,'XTick',[-3:9],'YTick',[0:0.1:1],'Yminortick','on');

text(yh.Extent(1),0.98*max(ylim),'C','Fontsize',32,'FontWeight','bold');

%% Plot temporal components
%% Set up information for the projection
%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
NatIData=NatIData(1:end-1);
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData'; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

% Load population density
load('PopulationDensity_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
P=log(P)';
% Load Rebel Control
load('RebelControl_Yemen.mat');
RC=RC';
% Load Health facility density
HS = shaperead([ pwd '\ShapeFile\healthsites.shp']); % Shape file for Yemen
load('PopulationSize_Yemen.mat');
H= GLevelHealthSites(HS,S);
H=10000.*H./AP';

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference



load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
NatIData=NatIData(1:end-1);
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

NW=length(NatIData); % number of week want to go out to
NWP=NW-length(WI(1,:));

% The size to separate the weeks in the x-label
    dW=4;
    % Set the X-tick labels to be Dates rather than numbers
    startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    %% Adjust ascpects of functions and data for the fitting

    %% Adjust ascpects of functions and data for the fitting

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
    %% Run the projection
    
    %% Run the logistic model with the data

    [Yt,~,~,~,~,~,~]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI));

    %% Run the projection
    
    temp=zeros(length(Yt(:,1)),1); % initialize the matrix for the projection
    Pt=zeros(length(Yt(:,1)),NWP); % used for the 
    for ii=1:NWP % Loop through the number of weeks that are to be projected
        WT=[WI(GNZI,:) temp]; % Need to append data to the end for the projection of incidence
        [temp2,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt,At,RCt]= LogisticModel(beta,WT,tA(GNZI,1:length(WT(1,:))),DB,DA,DAE,Ctv(GNZI,1:length(WT(1,:))),K,n,Rtv(GNZI,1:length(WT(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI)); % Run model with appendend data
        Pt(:,ii)=temp2(:,end); % Record the projection of incidence
        temp=[Pt(:,1:ii) zeros(length(Yt(:,1)),1)];  % temporary variable to be appended
    end
COther=(beta(2).*PDG+beta(3).*HFG+beta(4).*It+beta(9).*Gt+beta(11).*RCt);
    CRain=(beta(7).*IRt+beta(8).*Rt);
    CConflict=(beta(5).*IAt+beta(6).*ICt+beta(10).*At);
print(gcf,[pwd '\Figures\Figure2ABC.png'],'-dpng','-r600');
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.0488,0.575,0.9422,0.4]);
h=bar([(1+maxtau):152],[beta(8).*sum(Rt);beta(7).*sum(IRt)]','stacked','LineStyle','none'); 
h(2).FaceColor = 'flat';
    h(2).CData = hex2rgb('2171b5');
    h(1).FaceColor = 'flat';
    h(1).CData = hex2rgb('9ecae1');
box off;
ylabel({'Suspected cholera cases'});
legend({'Rainfall','Rainfall and Incidence'},'location','northwest');
legend boxoff;
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',10.^[0:5],'Fontsize',16,'XTickLabel',{''},'yscale','log');
ylim([1 100000]);
    xlim([1 NW]); % sets the x-limits of our x -axis
    xtickangle(45);
text(-5.536789293933579,82276.6933228582,'D','Fontsize',32,'FontWeight','bold');
subplot('Position',[0.0488,0.15,0.9422,0.4]);
h=bar([(1+maxtau):152],[beta(10).*sum(At);beta(5).*sum(IAt);beta(6).*sum(ICt)]','stacked','LineStyle','none'); 
h(3).FaceColor = 'flat';
    h(3).CData = hex2rgb('67000d');
h(2).FaceColor = 'flat';
    h(2).CData = hex2rgb('ef3b2c');
    h(1).FaceColor = 'flat';
    h(1).CData = hex2rgb('fc9272');
box off;
ylabel({'Suspected cholera cases'});
legend({'Attack','Attack and Incidence','Conflict and Incidence'},'location','northwest');
legend boxoff;
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',10.^[0:5],'Fontsize',16,'XTickLabel',XTL,'yscale','log');
ylim([1 100000]);
    xlim([1 NW]); % sets the x-limits of our x -axis
    xlabel('Week of report','Fontsize',24);
    xtickangle(45);
    text(-5.536789293933579,82276.6933228582,'E','Fontsize',32,'FontWeight','bold');
    print(gcf,[pwd '\Figures\Figure2DE.png'],'-dpng','-r600');
    
 %% Plot maps for rainfall and conflict/attacks
 
R=log10([sum(CRain,2)]+1);
CM=log10([sum(CConflict,2)]+1);
%MM=max([max(R) max(CM)]);
R=(R-min(R))./(max(R)-min(R));
CM=(CM-min(CM))./(max(CM)-min(CM));

N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});

Sm=cell(length(S),1); % cell array for the names of areas in shape file
for ii=1:length(Sm)
    Sm{ii}=S(ii).ADM1_EN; % record the names of the shape file
end


    f1=figure('units','normalized','outerposition',[0 0 1 1]);

    for ii=1:length(GNZI)    
        ff=find(contains(Sm,N(GNZI(ii)).G)); % see if the name is in the shape file based on the ordering OG
        mapshow(S(ff),'FaceColor',[0 0 0.6],'Edgecolor',[0 0 0],'LineWidth',1,'Facealpha',R(ii));hold on % shape the area in with color C(ii,:)
    end

    mapshow(S(12),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)

    mapshow(S(21),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
    xl=xlim; % saves the xlim to xl
    yl=ylim;% saves the ylim to yl

    % Produces a cross hatch for polygon 12
    CrossHatchMap(min(S(12).X),max(S(12).X),min(S(12).Y),max(S(12).Y),0.01,[0 0 0],S(12))
    % Produces a cross hatch for polygon 21
    CrossHatchMap(min(S(21).X),max(S(21).X),min(S(21).Y),max(S(21).Y),0.025,[0 0 0],S(21))
    xlim(xl); % set xlimits to xl
    ylim(yl); % set ylimits to yl
    set(gca,'visible','off') % remove the axis from the map 
    
        text(min(xl),max(yl),'F','Fontsize',32,'Fontweight','bold');
      
    print(gcf,[pwd '\Figures\Figure2F.png'],'-dpng','-r600');

% 
% 
% %% Govnorate level various waves and the impact of conflict
% 
% 
% 
    f1=figure('units','normalized','outerposition',[0 0 1 1]);
    for ii=1:length(GNZI)    
        ff=find(contains(Sm,N(GNZI(ii)).G)); % see if the name is in the shape file based on the ordering OG
        mapshow(S(ff),'FaceColor',[0.6 0 0],'Edgecolor',[0 0 0],'LineWidth',1,'Facealpha',CM(ii));hold on % shape the area in with color C(ii,:)
    end

    mapshow(S(12),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)

    mapshow(S(21),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
    xl=xlim; % saves the xlim to xl
    yl=ylim;% saves the ylim to yl

    % Produces a cross hatch for polygon 12
    CrossHatchMap(min(S(12).X),max(S(12).X),min(S(12).Y),max(S(12).Y),0.01,[0 0 0],S(12))
    % Produces a cross hatch for polygon 21
    CrossHatchMap(min(S(21).X),max(S(21).X),min(S(21).Y),max(S(21).Y),0.025,[0 0 0],S(21))
    xlim(xl); % set xlimits to xl
    ylim(yl); % set ylimits to yl
    set(gca,'visible','off') % remove the axis from the map 
    
        text(min(xl),max(yl),'G','Fontsize',32,'Fontweight','bold');        
     print(gcf,[pwd '\Figures\Figure2G.png'],'-dpng','-r600');
