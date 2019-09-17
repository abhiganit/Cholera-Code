
load('Yemen_Gov_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

WI=IData'; % Transpose the data set such that the number of areas is the row
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection

AF=2;
XU=ones(1,11);
load([pwd '\Tables\TestProjectionModelGA-XU=2047-CF=2-RIF=0-PF=1-tau=1  1  1  2  1  4  1  1  3  1.mat'],'par');


%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
    C=ImpactConflict([0:10],K,n,CF);
subplot(3,2,1);bar([0:10],C)
A=ImpactAttack([zeros(1,21) 1 zeros(1,21)],0,DAE,0,0);
subplot(3,2,2);bar([-21:21],A); 

A=ImpactAttack([zeros(1,21) 1 zeros(1,21)],DB,DA,0,0);
subplot(3,2,3);bar([-21:21],A); 

R=ImpactRainfall(linspace(0,16,10001),RIF,rl);
subplot(3,2,4);plot(linspace(0,16,10001),R); 


R=ImpactRainfall(linspace(0,16,10001),RF,rh);
subplot(3,2,5);plot(linspace(0,16,10001),R); 