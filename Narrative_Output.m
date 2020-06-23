[WI2,~,tA2,Rtv2,~,P2,RC2,H2,WPINm2,FPINm2,Dieselt2,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
[WIG2,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,CIG] = LoadYemenData; % Load the data used to construct the figure
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

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

% Record the names for the IDP calculation
Sm=cell(length(S),1); % allocate space
for ii=1:length(Sm)
Sm{ii}=S(ii).ADM1_EN; % record name
end

SD = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen

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

load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,SD,153);

temp=GLevelConflict(ProC,S([9 2]),153);
Ctv=[Ctv; sum(Ctv(1:3,:),1);temp];

load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,SD,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
temp=GLevelConflict(YASt,S([9 2]),153);
Mt=[Mt; sum(Mt(1:3,:),1);temp];


IndW=[1 21; 22 74; 75 121; 122 149]+maxtau; % Index of wave for the data used in the regression model
IndW(1,1)=1;

IW=7.*(([1; 22 ; 75 ; 122; 150]-1)+maxtau); % The 150 is the start of the week we do not have data for and we are subtracting a week for the index of the week as the index zero is Oct 3, 2016
IW=[IW(1) IW(2)-1 IW(2) IW(3)-1 IW(3) IW(4)-1 IW(4) IW(5)-1];
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+IW],'mmm.dd,yyyy');
%% Amanat Al Asimah

fprintf(['Number of suspected cases in Amanat Al Asimah from ' datestr(startDateofSim) ' to ' datestr(startDateofSim+7*(IndW(2,1)-1)) ': %d \n'], CIG(9,IndW(1,2)))
fprintf(['Number of conflict evetns in Amanat Al Asimah from ' datestr(startDateofSim) ' to ' datestr(startDateofSim+7*(IndW(2,1)-1)) ': %d \n'], sum(Ctv(end-1,IndW(1,1):IndW(1,2))))
fprintf(['Number of shelling and air attacks in Amanat Al Asimah from ' datestr(startDateofSim) ' to ' datestr(startDateofSim+7*(IndW(2,1)-1)) ': %d \n'], sum(Mt(end-1,IndW(1,1):IndW(1,2))))
