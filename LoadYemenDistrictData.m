function [WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,VT1,VT2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData
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
% maxtau - the maximum lag considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

load('Yemen_Gov_Incidence.mat'); % Incidence data
CItemp=cumsum(IData',2);
CItemp=[zeros(size(CItemp(:,1))) CItemp(:,1:end-1)];
CItemp=CItemp(:,TruncV:end);
load('Yemen_District_Incidence.mat'); % Incidence data
CI=cumsum(IData',2);
CI=[zeros(size(CI(:,1))) CI(:,1:end-1)];
CI(end-1,:)=CItemp(9,:);
CI(end,:)=CItemp(2,:);
S = shaperead([ pwd '/ShapeFile/yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

% Record the names for the IDP calculation
Sm=cell(length(S),1); % allocate space
for ii=1:length(Sm)
Sm{ii}=S(ii).ADM1_EN; % record name
end

SD = shaperead([ pwd '/ShapeFile/yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen

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


SD=SD([29 31 71 fS' fA']);
% Record the names for the IDP calculation
SDm=cell(length(SD),1); % allocate space
SDGm=cell(length(SD),1); % allocate space
for ii=1:length(SDm)
SDm{ii}=SD(ii).ADM2_EN; % record name
SDGm{ii}=SD(ii).ADM1_EN; % record name
end

WI=IData'; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,SD,153);

temp=GLevelConflict(ProC,S([9 2]),153);
Ctv=[Ctv; sum(Ctv(1:3,:),1);temp];
Ctv=Ctv(:,TruncV:end); % Incidence data for the districts starts at may 1 2017
load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,SD,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
temp=GLevelConflict(YASt,S([9 2]),153);
Mt=[Mt; sum(Mt(1:3,:),1);temp];
Mt=Mt(:,TruncV:end); % Incidence data for the districts starts at may 1 2017
% Use attacks with "water in description"
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,SD,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
temp=GLevelConflict(ProA,S([9 2]),153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
tA=[tA; sum(tA(1:3,:),1);temp];
tA=tA(:,TruncV:end); % Incidence data for the districts starts at may 1 2017
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Ginx=zeros(length(SDm)+3,1);
for ii=1:length(SDm)  
   Ginx(ii)=find(strcmp(SDGm(ii),Sm));
end
Ginx(length(SDm)+1)=Ginx(1);
Ginx(length(SDm)+2)=9;
Ginx(length(SDm)+3)=2;
Rtv=Rtv(Ginx,TruncV:153); % truncate the railfall data to the time period spceified
% W = shaperead([ pwd '/ShapeFile/Wadies.shp']); % Shape file for Yemen water course
% Mt=WaterCourseConnection(W,S); % Calculate the contact matrix for the water course
% Load population density
load('Area_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
Atemp=A(Ginx);
for ii=1:length(SDm)    
   Atemp(ii)=Atemp(ii).*SD(ii).Shape_Area./(S(Ginx(ii)).Shape_Area);
end
Atemp(length(SDm)+1)=sum(Atemp(1:3));
A=Atemp;
load('PopulationSize_DistrictYemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
P=log([ repmat(AP(:,1)./A,1,NW2016) repmat(AP(:,2)./A,1,52)  repmat(AP(:,3)./A,1,52)  repmat(AP(:,4)./A,1,NW2019)]);
P=P(:,TruncV:end); % Incidence data for the districts starts at may 1 2017
% Load Rebel Control
load('RebelControl_Yemen.mat');
RC=RC';
RC=RC(Ginx);
% Load Health facility density
HS = shaperead([ pwd '/ShapeFile/healthsites.shp']); % Shape file for Yemen
load('PopulationSize_DistrictYemen.mat');
H= GLevelHealthSites(HS,SD);
temp=GLevelHealthSites(HS,S([9 2]));
H=[H; sum(H([1:3]));temp];
H=10000.*[ repmat(H./AP(:,1),1,NW2016) repmat(H./AP(:,2),1,52)  repmat(H./AP(:,3),1,52)  repmat(H./AP(:,4),1,NW2019)];
H=H(:,TruncV:end); % Incidence data for the districts starts at may 1 2017
% WASH People in Need Desnity
load('WASH_PIN_DistrictYemen.mat');
WPINm=[ repmat(WPIN(:,1)./AP(:,1),1,NW2016) repmat(WPIN(:,2)./AP(:,2),1,52)  repmat(WPIN(:,3)./AP(:,3),1,52)  repmat(WPIN(:,4)./AP(:,4),1,NW2019)];
WPINm=WPINm(:,TruncV:end); % Incidence data for the districts starts at may 1 2017

load('Food_PIN_DistrictYemen.mat');
FPINm=[ repmat(FPIN(:,1)./AP(:,1),1,NW2016) repmat(FPIN(:,2)./AP(:,2),1,52)  repmat(FPIN(:,3)./AP(:,3),1,52)  repmat(FPIN(:,4)./AP(:,4),1,NW2019)];
FPINm=FPINm(:,TruncV:end); % Incidence data for the districts starts at may 1 2017

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end); % Incidence data for the districts starts at may 1 2017


%% Diesel and Wheat prices

load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Diesel_Gov_Yemen.mat')
load('Wheat_Gov_Yemen.mat')
Diesel=Diesel'-min(Diesel(Diesel>0));
Wheat=Wheat'-min(Wheat(Wheat>0));
Wheatt=zeros(size(IData'));
Dieselt=zeros(size(IData'));

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


% Temprature
Temptv=zeros(size(WI));
load('Temprature_Gov.mat','Temp_G');
EndFirstEpiWeek= datenum('10-09-2016');% Start of second epiweek is oct 10, 2016
for ii=1:length(IData(:,1))
    f=find(EndFirstEpiWeek+7*(ii-1)<=endData,1);% Need the first one to satisfy this condition as we increase over time and will elminate the other past months
    Wheatt(:,ii)=Wheat(:,f);
    Dieselt(:,ii)=Diesel(:,f);
    Temptv(:,ii)=Temp_G(:,month(EndFirstEpiWeek+7*(ii-1)));
end

Wheattemp=Wheatt(:,TruncV:end);
Dieseltemp=Dieselt(:,TruncV:end);
Temptvtemp=Temptv(:,TruncV:end);

Wheatt=Wheattemp([5 5 5 9.*ones(1,length(fS)) 2.*ones(1,length(fA)) 5 9 2],:);
Dieselt=Dieseltemp([5 5 5 9.*ones(1,length(fS)) 2.*ones(1,length(fA)) 5 9 2],:);
Temptv=Temptvtemp([5 5 5 9.*ones(1,length(fS)) 2.*ones(1,length(fA)) 5 9 2],:);
%% Comput the attack rate per 10,000
WI=10000.*WI./PopS;

Ctv=log(1+Ctv./repmat(A,1,length(Ctv(1,:))).*PopS);
Mt=log(1+Mt./repmat(A,1,length(Mt(1,:))).*PopS);

% Vaccination doses 
[V1,V2]= VaccinationTime(2,153);

V1=V1(:,TruncV:end);
V2=V2(:,TruncV:end);

VT1=V1./PopS;
VT2=V2./PopS;

GV=sum(VT1+VT2,2);
GV(GV>0)=1;

%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model


end

