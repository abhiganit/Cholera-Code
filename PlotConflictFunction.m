
load('Yemen_Gov_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

WI=IData'; % Transpose the data set such that the number of areas is the row
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection

AF=2;
XU=ones(1,11);
load([pwd '\Tables\ProjectionModelGA-XU=2047-CF=1-RIF=1-PF=1-tau=1  1  1  3  3  4  1  1  2  1.mat'],'par');


%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterGA(par,XU);
    
C=ImpactConflict(Ctv(GNZI,(1+maxtau-tau(5)):(end-tau(5))),K,n,CF);
subplot(3,1,1);plot((1+maxtau):(131),sum(C))
ylabel('Conflict function');
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,131);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
A=ImpactAttack(tA(GNZI,:),0,DAE,tau(9),maxtau);
AF=[zeros(1,maxtau) sum(A)];
ST=sum(tA(GNZI,:));
A=ImpactAttack(tA(GNZI,:),DB,DA,tau(4),maxtau);

AF2=[zeros(1,maxtau) sum(A)];
ST=sum(tA);
subplot(3,1,2);bar([1:131],[ST]); 
ylabel('Number of attacks');

subplot(3,1,3);bar([1:131],[AF;AF2]','stacked'); 
ylabel('Attack functions');
legend({'Attack','Attack and Incidence'})