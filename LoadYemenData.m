function [WI,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,VT1,VT2,GNZI,GV,maxtau] = LoadYemenData
% Loads the Data needed to rum the regression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WI - Weekly incidence (22x153)
% Ctv - Weekly level of conflict (22x153)
% tA - Weekly number of attacks (22x153)
% Rtv - Weekly level of rainfall (22x153)
% Mt - Water course connection (22x22)
% P - log of population density (22x1)
% RC - Rebel control (22x1)
% H - Health sites per 10000 (22x1)
% WPIN - transformation of density of WASH people in need (22x1)
% Et - External incidence due to IDP 
% GNZI - Gov with non-zero incidence
% GV - Inidciator if vaccination occurred in governorate
% maxtau - the maximum lag considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Yemen_Gov_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

% Record the names for the IDP calculation
Sm=cell(length(S),1); % allocate space
for ii=1:length(Sm)
Sm{ii}=S(ii).ADM1_EN; % record name
end

WI=IData'; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
% Use attacks with "water in description"
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified
% W = shaperead([ pwd '\ShapeFile\Wadies.shp']); % Shape file for Yemen water course
% Mt=WaterCourseConnection(W,S); % Calculate the contact matrix for the water course
% Load population density
load('Area_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
load('PopulationSize_Yemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
P=log([ repmat(AP(:,1)./A,1,NW2016) repmat(AP(:,2)./A,1,52)  repmat(AP(:,3)./A,1,52)  repmat(AP(:,4)./A,1,NW2019)]);
% Load Rebel Control
load('RebelControl_Yemen.mat');
RC=RC';
% Load Health facility density
HS = shaperead([ pwd '\ShapeFile\healthsites.shp']); % Shape file for Yemen
load('PopulationSize_Yemen.mat');
H= GLevelHealthSites(HS,S);
H=10000.*[ repmat(H./AP(:,1),1,NW2016) repmat(H./AP(:,2),1,52)  repmat(H./AP(:,3),1,52)  repmat(H./AP(:,4),1,NW2019)];
% WASH People in Need Desnity
load('WASH_PIN_Yemen.mat');
WPINm=[ repmat(WPIN(:,1)./AP(:,1),1,NW2016) repmat(WPIN(:,2)./AP(:,2),1,52)  repmat(WPIN(:,3)./AP(:,3),1,52)  repmat(WPIN(:,4)./AP(:,4),1,NW2019)];


load('Food_PIN_Yemen.mat');
FPINm=[ repmat(FPIN(:,1)./AP(:,1),1,NW2016) repmat(FPIN(:,2)./AP(:,2),1,52)  repmat(FPIN(:,3)./AP(:,3),1,52)  repmat(FPIN(:,4)./AP(:,4),1,NW2019)];
%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
% load('Yemen_IDP.mat'); % load IDP cell matrix
% IDPtAbs=IDPMatrix(IDP,Sm,153); % Calculate temporal IDP absolute numbers
% IDPt = Et(WI,PopS,IDPtAbs);
%% Comput the attack rate per 10,000
WI=10000.*WI./PopS;

%% Diesel and Wheat prices
load('Diesel_Gov_Yemen.mat')
load('Wheat_Gov_Yemen.mat')
Diesel=Diesel';
Wheat=Wheat';
Wheatt=zeros(size(WI));
Dieselt=zeros(size(WI));

% The dates for the prices being used
startData = [ ];
endData = [ ];
for yy=2016:2019
    if(yy==2016)
        for mm=10:12
            startData=[startData; datenum(yy,mm,1)];
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end
        end
    elseif(yy==2019)        
        for mm=1:9
            startData=[startData; datenum(yy,mm,1)];    
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end        
        end
    else
        for mm=1:12
            startData=[startData; datenum(yy,mm,1)];
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end
        end
    end
end

EndFirstEpiWeek= datenum('10-09-2016');% Start of second epiweek is oct 10, 2016
for ii=1:length(WI(1,:))
    f=find(EndFirstEpiWeek+7*(ii-1)<=endData,1);% Need the first one to satisfy this condition as we increase over time and will elminate the other past months
    Wheatt(:,ii)=Wheat(:,f);
    Dieselt(:,ii)=Diesel(:,f);
end

% Vaccination doses 
[V1,V2]= VaccinationTime(1,153);

VT1=V1./PopS;
VT2=V2./PopS;

GV=sum(VT1+VT2,2);
GV(GV>0)=1;

%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model


end

