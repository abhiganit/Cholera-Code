%% Propogates from the start of the epidemic

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

maxtau=4;

load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,131);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

NW=length(NatIData); % number of week want to go out to
NWP=131-maxtau;

CInc=zeros(20,126);
CRain=zeros(20,126);
CConflict=zeros(20,126);
OvC=zeros(1,11);
    AF=2;
    XU=ones(1,11);%[1 1 1 1 0 0 1 1 1 0 0];%ones(1,11);
    load([pwd '\Tables\BackSelectModel-PercentDataSet=90-XU=2047-CF=2-RIF=1-PF=2-tau=1  1  1  2  1  2  1  1  3  1.mat'],'par');
    %% Adjust ascpects of functions and data for the fitting

    %% Adjust ascpects of functions and data for the fitting

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
    %% Run the projection
    
    temp=zeros(length(GNZI),1); % initialize the matrix for the projection
    Pt=zeros(length(GNZI),NWP); % used for the 
    for ii=1:NWP % Loop through the number of weeks that are to be projected
        WT=[WI(GNZI,1:maxtau) temp]; % Need to append data to the end for the projection of incidence
        [temp2,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt,At,RCt]= LogisticModel(beta,WT,tA(GNZI,1:length(WT(1,:))),DB,DA,DAE,Ctv(GNZI,1:length(WT(1,:))),K,n,Rtv(GNZI,1:length(WT(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI)); % Run model with appendend data
        Pt(:,ii)=temp2(:,end); % Record the projection of incidence
        temp=[Pt(:,1:ii) zeros(length(GNZI),1)];  % temporary variable to be appended
    end
    plot([1:NWP]+maxtau,sum(Pt),'k','LineWidth',2); hold on;
    
    XU([5 6 10 11])=0;
    temp=zeros(length(GNZI),1); % initialize the matrix for the projection
    Pt=zeros(length(GNZI),NWP); % used for the 
    for ii=1:NWP % Loop through the number of weeks that are to be projected
        WT=[WI(GNZI,1:maxtau) temp]; % Need to append data to the end for the projection of incidence
        [temp2,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt,At,RCt]= LogisticModel(beta.*XU,WT,tA(GNZI,1:length(WT(1,:))),DB,DA,DAE,Ctv(GNZI,1:length(WT(1,:))),K,n,Rtv(GNZI,1:length(WT(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI)); % Run model with appendend data
        Pt(:,ii)=temp2(:,end); % Record the projection of incidence
        temp=[Pt(:,1:ii) zeros(length(GNZI),1)];  % temporary variable to be appended
    end
    plot([1:NWP]+maxtau,sum(Pt),'b','LineWidth',2);
    Pt2=sum(Pt);
    XU([5 6 10 11])=1;
    XU([ 7 8])=0;
    temp=zeros(length(GNZI),1); % initialize the matrix for the projection
    Pt=zeros(length(GNZI),NWP); % used for the 
    for ii=1:NWP % Loop through the number of weeks that are to be projected
        WT=[WI(GNZI,1:maxtau) temp]; % Need to append data to the end for the projection of incidence
        [temp2,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt,At,RCt]= LogisticModel(beta.*XU,WT,tA(GNZI,1:length(WT(1,:))),DB,DA,DAE,Ctv(GNZI,1:length(WT(1,:))),K,n,Rtv(GNZI,1:length(WT(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI)); % Run model with appendend data
        Pt(:,ii)=temp2(:,end); % Record the projection of incidence
        temp=[Pt(:,1:ii) zeros(length(GNZI),1)];  % temporary variable to be appended
    end
    plot([1:NWP]+maxtau,sum(Pt),'r','LineWidth',2);

    plot([1:NWP]+maxtau,sum(Pt)+Pt2,'color',[0.6 0.6 0.6],'LineWidth',2);
