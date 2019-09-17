close all;
% F=struct('N',{'None','Low','High','LowHigh'});
% CFF=struct('N',{'None','Linear','Linear_Threshold','Poly_Threshold'});
% AFF=struct('N',{'None','Before','After','BeforeAfter'});
% 
% T=readtable([pwd '\Tables\ProjectionModelFitSummary.dat']); 
% XT=[T.Intercept    T.Population    T.HealthFacilities    T.Incidence    T.Attack_Incidence    T.Conflcit_Incidence    T.Rainfall_Incidence    T.Precipitation    T.External_Incidence    T.Attack T.RebelControl];
% AIC=T.AIC;
% dAIC=AIC-min(AIC);
% w=exp(-dAIC./2)./sum(exp(-dAIC./2));
% NN=length(w);
%% Set up information for the projection
%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData'; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,131);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
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
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,131);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

NW=length(NatIData); % number of week want to go out to
NWP=NW-length(WI(1,:));

CInc=zeros(20,126);
CRain=zeros(20,126);
CConflict=zeros(20,126);
OvC=zeros(1,11);
% Run best model
% indx=find(dAIC==0);
% indx=96;% indx(1);
%     XU=XT(indx,:);
%     tau=XU(2:end);
%     XU(XU>1)=1;
%     AF=-2;
%     ftemp=[];
%     while (isempty(ftemp))
%         AF=AF+1;
%         ftemp=find(contains(T.Attack_Function(indx),AFF(AF+2).N));
%     end
%     CF=-2;
%     ftemp=[];
%     while (isempty(ftemp))
%         CF=CF+1;
%         ftemp=find(contains(T.Conflict_Function(indx),CFF(CF+2).N));
%     end
%     
%     RIF=-2;
%     ftemp=[];
%     while (isempty(ftemp))
%         RIF=RIF+1;
%         ftemp=find(contains(T.Rainfall_Incidence_Function(indx),F(RIF+2).N));
%     end
%     
%     RF=-2;
%     ftemp=[];
%     while (isempty(ftemp))
%         RF=RF+1;
%         ftemp=find(contains(T.Precipitation_Function(indx),F(RF+2).N));
%     end
%     
    AF=2;
    XU=ones(1,11);
    load([pwd '\Tables\TestProjectionModelGA-XU=2047-CF=2-RIF=0-PF=1-tau=1  1  1  2  1  4  1  1  3  1.mat'],'par');
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

    %Plot national data in a scatter plot
    scatter([1:NW],NatIData,40,'k','filled'); hold on;
    %plot the model fit
    M=sum(Yt,1); % Aggregate the incidence in the various areas for the fitted portion
    NWF=length(WI(1,:)); % Number of weeks for the fitting
     plot([(1+maxtau):NWF],M,'b','LineWidth',2);  % Plot the fitted portion of the mode

     %Plot the model projection
    M=sum([Yt(:,end) Pt],1); % Aggregate the incidence in the various areas for the projection
    plot([0:length(M)-1]+NWF,M,'r','LineWidth',2); hold off;% Plot the projected portion of the mode

    % Adjust characteristics of the figure
     box off; % removes the outside box on the figure
    xlim([1 NW]); % sets the x-limits of our x -axis
    ylim([0 55000]); %sets the y-limts of the y-axis
    % The size to separate the weeks in the x-label
    dW=4;
    % Set the X-tick labels to be Dates rather than numbers
    startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    % changing the aspects of the axis for the the current figure 
    set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:5000:55000],'Fontsize',16,'XTickLabel',XTL);
    % Rotates the xticklabel 
    xtickangle(45) 
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    % Xlable of the figure
    xlabel('Week of report','Fontsize',24);
    % ylable of the figure
    ylabel('Number of reported cases','Fontsize',24);
    % Puts text in the figure for labelling the fit and projection
    text(NW-14,52500,'Model fit','color','b','Fontsize',18); 
    text(NW-14,50000,'Model projection','color','r','Fontsize',18);
figure('units','normalized','outerposition',[0 0 1 1]);

   subplot(3,1,1);plot((1+maxtau):(130),beta(6).*sum(ICt))
ylabel('Suspected cholera cases');
ST=sum(tA(GNZI,:));
subplot(3,1,2);bar([1:131],[ST]); 
ylabel('Number of attacks');

subplot(3,1,3);bar([(1+maxtau):130],[beta(10).*sum(At);beta(5).*sum(IAt)]','stacked'); 
ylabel('Suspected cholera cases');
legend({'Attack','Attack and Incidence'})