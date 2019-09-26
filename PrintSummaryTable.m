%% Prints a summary of the information for the different govnorates

load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

% Load population density
load('PopulationDensity_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
%P=log(P)';
% Load Rebel Control
load('RebelControl_Yemen.mat');
RC=RC';
% Load Health facility density
HS = shaperead([ pwd '\ShapeFile\healthsites.shp']); % Shape file for Yemen
load('PopulationSize_Yemen.mat');
H= GLevelHealthSites(HS,S);
H=10000.*H./AP';
RCYN=struct('B',{'No','Yes'});

f1=fopen('LatexTableGov.txt','w');

fprintf(f1,'\\begin{table} \n \\caption{A summary of the govnorate level information used in the analysis of the cholera outbreak in Yemen} \n \\begin{tabular}{cccccccc} \n');
fprintf(f1,'Govnorate & Population density (per km$^2$) & Health sites per 10,000 & Total suspected cases & Avg. Weekly Rainfall (mm) & Avg. Weekly number of conflict & Total number of attacks & Under rebel control \\\\ \n \\hline ');
for ii=1:length(S)-1
    fprintf(f1,[S(ii).ADM1_EN ' & %3.2E & %3.2f & %3.2E & %4.2f & %4.2f & %d & ' RCYN(RC(ii)+1).B ' \\\\ \n'],[P(ii) H(ii) sum(WI(:,ii)) mean(Rtv(ii,:)) mean(Ctv(ii,:)) sum(tA(ii,:)) ]);
end
ii=length(S);
fprintf(f1,[S(ii).ADM1_EN ' & %3.2E & %3.2f & %3.2E & %4.2f & %4.2f & %d & ' RCYN(RC(ii)+1).B ' \n'],[P(ii) H(ii) sum(WI(:,ii)) mean(Rtv(ii,:)) mean(Ctv(ii,:)) sum(tA(ii,:))  ]);
fprintf(f1,'\\end{tabular} \n \\end{table}');
fclose(f1);
