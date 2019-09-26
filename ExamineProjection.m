close all;
%% Set up information for the projection
%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
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
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

NW=length(NatIData)-1; % number of week want to go out to
NWP=NW-length(WI(1,:));

    AF=2;
    XU=[1 1 1 1 0 0 1 1 1 0 0];%ones(1,11);%
    load([pwd '\Tables\BackSelectModel-PercentDataSet=100-XU=463-CF=0-RIF=0-PF=0-tau=1  1  1  4  4  1.mat'],'par');
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

figure('units','normalized','outerposition',[0 0 1 1]);

    subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
   PDS=0.9;
    %Plot national data in a scatter plot
    scatter([1:floor(length(NatIData)*PDS)],NatIData(1:floor(length(NatIData)*PDS)),40,'k','filled'); hold on;
    
    scatter([floor(length(NatIData)*PDS)+1:NW],NatIData([floor(length(NatIData)*PDS)+1:NW]),40,'k','filled'); hold on;
    %plot the model fit
    Mt=sum([Yt Pt]);
     plot([(1+maxtau):floor(length(NatIData)*PDS)],Mt([(1+maxtau):floor(length(NatIData)*PDS)]-maxtau),'k','LineWidth',2);  % Plot the fitted portion of the mode

     %Plot the model projection
    plot([floor(length(NatIData)*PDS):NW],Mt(floor(length(NatIData)*PDS)-maxtau:end),'k','LineWidth',2); hold off;% Plot the projected portion of the mode

    % Adjust characteristics of the figure
     box off; % removes the outside box on the figure
    xlim([1 NW]); % sets the x-limits of our x -axis
    ylim([0 65000]); %sets the y-limts of the y-axis
    % The size to separate the weeks in the x-label
    dW=4;
    % Set the X-tick labels to be Dates rather than numbers
    startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    % changing the aspects of the axis for the the current figure 
    set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:5000:65000],'Fontsize',16,'XTickLabel',XTL);
    % Rotates the xticklabel 
    xtickangle(45) 
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    % Xlable of the figure
    %xlabel('Week of report','Fontsize',24);
    % ylable of the figure
    ylabel('Number of suspected cholera cases','Fontsize',24);
    % Puts text in the figure for labelling the fit and projection
%     text(NW-14,52500,'Model fit','color','b','Fontsize',18); 
%     text(NW-14,50000,'Model projection','color','r','Fontsize',18); 
hold on;

 AF=2;
    XU=ones(1,11);%
    load([pwd '\Tables\BackSelectModel-PercentDataSet=100-XU=2047-CF=1-RIF=0-PF=1-tau=1  1  1  3  1  4  4  1  2  1.mat'],'par');
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
    
    PDS=0.9;
    %plot the model fit
    Mt=sum([Yt Pt]);
     plot([(1+maxtau):floor(length(NatIData)*PDS)],Mt([(1+maxtau):floor(length(NatIData)*PDS)]-maxtau),'r','LineWidth',2);  % Plot the fitted portion of the mode

     %Plot the model projection
    plot([floor(length(NatIData)*PDS):NW],Mt(floor(length(NatIData)*PDS)-maxtau:end),'r','LineWidth',2); hold off;% Plot the projected portion of the mode

figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.0708,0.5901,0.897162184873949,0.4]);
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

subplot('Position',[0.0708,0.166,0.897162184873949,0.4]);
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