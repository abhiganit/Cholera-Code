% Close all figures and clear all other data and clear command window
close all;
clear;
clc;
%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData'; % Transpose the data set such that the number of areas is the row
Ctv=GLevelConflict(ProC,S,length(WI(1,:))); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
%% NEED TO ENSURE THAT THE ROWS ARE THE PROPER GOVERNORATES %%
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas
Rtv=Rtv(:,1:length(WI(1,:))); % truncate the railfall data to the time period spceified
load('Attack_Yemen_Time_Location_Project_Forward.mat');
tA=GLevelConflict(ProA,S,length(WI(1,:)));% % Send attack data (time, latitude, longatude) and shapefile of area wanted to catagorize
tA(tA>1)=1; % Find the weeks where there are more than one attack that occurs and assume that there is only a single attack

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

%% Specify ascpects of functions and data for the fitting
tau=[1 1 1 1 2]; % Specify a lag of all the factors that are integrated in the model

%Specify what you want to use in in the regression model

XU=[0 1 1 1 0 1];
f=find(XU(2:end)==0);
tau(f)=0; % set the lag for components not included to zero

% Specify the environmental factors that you want to use in the model fitting
A=1; % A=1 the effects of water related attacks will be used (A=0 not used)
C=1; % C=1 the effects of water related attacks will be used (C=0 not used)
R=0; % R=1 the effects of water related attacks will be used (R=0 not used)

%Specify the attack function to be used
AF=1; % AF=0 attack only has effect before; AF=1 Attack has effect only after; AF=2; Attack has effect before and after

%Specify the conflict function to be used
CF=2; % CF=0 linear effect; CF=1 Hill function with n=1; CF=2; Full hill function

% Specify the rainfall function to be used
RF=2; % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall

% Specify the lower bounds for the estimated parameters
% beta(1:6)=XU.*[x(1) 10.^(x(2:6))]
% 
% %Attack associated paramters
% DB=10.^x(7);
% DA=10.^x(8);
% % Conflict associated paramters
% K=10.^x(9);
% n=10.^x(10);
% % Rainfall assocaited paramters
% rl=10.^x(11);
% rh=10.^x(12);
% Environmental
% ma=10.^x(13);
% mc=10.^x(14);
% mr=10.^x(15);


lb=[-inf.*ones(1,12) -3 -3 -3]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=log10([inf inf inf inf inf inf 1 1 inf inf 10 10 1 1 1]); % specify the upperbound for the parameters 


%% Run the fitting algorithm
NS=10000; % The number of initial starting points used in the fitting process
opts= optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-12),'TolX',10^(-12),'UseParallel',true,'Display','off');  % Specifies the conditions that will be used in the fitting processs
problem = createOptimProblem('lsqnonlin','objective',@(x)OFunc(x,WI(GNZI,:),A,C,R,tA(GNZI,:),Ctv(GNZI,:),Rtv(GNZI,:),XU,tau,AF,CF,RF),'x0',log10(rand(15,1)),'lb',lb,'ub',ub,'Aineq',[],'bineq',[],'Aeq',[],'beq',[],'options',opts);
ms = MultiStart('UseParallel','always'); % specifies that we run the algorithm in parallel
[par,fval] = run(ms,problem,NS); %starts at NS random initial points to thoroghly search the paramter space
% Evaluate the number of paramters that are being used in the estimation 
[k,beta,DB,DA,K,n,rl,rh,ma,mc,mr]=RetParameter(par,XU,A,C,R,AF,CF,RF);

% Computing the statistical of the bets estimates for the regresison model
    f=find(XU==1); %  fiind components that are included
    be=beta(f); % increase index by one as we do not remove the inital beta_0
    [Yt,It,IAt,ICt,IRt,IEt]=LogisticModel(beta,WI(GNZI,:),A,tA(GNZI,:),DB,DA,C,Ctv(GNZI,:),K,n,R,Rtv(GNZI,:),RF,rl,rh,tau,CF,ma,mc,mr); % Returns the incidence in matrix form of size Ng X (NW-tau)
    resid=OFunc(par,WI,A,C,R,tA,Ctv,Rtv,XU,tau,AF,CF,RF); % Returns the vector of residulas for each of the data points
    X=[]; % Construct the input vector based on the variables used in the fitting
    twotail=[]; % for the two-tail or single tail p-value (1 - uses two-tail)
    % If the constant is used
    if(XU(1)==1)
        X=[X;ones(1,length(It(:)))]; % add the constant one to the input
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    %If past incidence is used
    if(XU(2)==1)
        X=[X It(:)]; % add the past incidence
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    %If attacks are used
    if(XU(3)==1)
        X=[X IAt(:) ]; % add product of past incidence and attacks
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If confclit is used
    if(XU(4)==1)
        X=[X ICt(:)]; % add product of past incidence and conflict
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If rainfall is used
    if(XU(5)==1)
        X=[X IRt(:)]; % add product of past incidence and rainfall
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If environmental factors are used
    if(XU(6)==1)
        X=[X IEt(:)]; % add product of past incidence and environmental factors
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % Calculate the statistics for for the coefficients of the model
    [SE,tStat,pValue] = LinRegressionStat(be,resid,k,X,twotail);   % We have specified the two-tail p-value as seen above
    
    % Constrcuts a table that is saved an can be viewed later
    f1 = fopen(['Table-XU=' num2str(XU(1)+2.*XU(2)+4.*XU(3)+8.*XU(4)+16.*XU(5)+32.*XU(6)) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-C=' num2str(C) '-A=' num2str(A) '-R=' num2str(R) '.dat'],'w'); % Writing to the file name
    CN=struct('N',{'Intercept','Incidence','Attack','Conflcit','Rainfall','Environment'}); % The name of the variables that we can call
    fprintf(f1, 'Coefficient \t Value \t Stand._Error \t t-Stat \t p-value \n'); % The header for the table
    for ii=1:6 % Go through all variables to see which are included
        if(XU(ii)~=0) % print on the variable that are included in the model
            fprintf(f1, [CN(ii).N ' \t %4.3f \t %4.3f \t %4.3f \t %3.2E \n'],[be(sum(XU(1:ii))) SE(sum(XU(1:ii))) tStat(sum(XU(1:ii))) pValue(sum(XU(1:ii)))]); % writes the information to the file
        end
    end
    fclose(f1); % closes the file
    readtable(['Table-XU=' num2str(XU(1)+2.*XU(2)+4.*XU(3)+8.*XU(4)+16.*XU(5)+32.*XU(6)) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-C=' num2str(C) '-A=' num2str(A) '-R=' num2str(R) '.dat']) % reads and prints the table to the command window
    save(['Model-XU=' num2str(XU(1)+2.*XU(2)+4.*XU(3)+8.*XU(4)+16.*XU(5)+32.*XU(6)) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-C=' num2str(C) '-A=' num2str(A) '-R=' num2str(R) '.mat'],'k','beta','DB','DA','K','n','rl','rh','ma','mc','mr','fval','tau'); % Saves the information for the specified model for can load and run model after

%% Plot the enviromental functions for the model


figure('units','normalized','outerposition',[0 0 1 1]);

% Dimensions for the subpanels 
hei=0.36; % heigh of figure
wid=0.42; % width of figure
xxs=linspace(0.075,0.99-0.42,2); % x-location
yys=linspace(0.97-hei,0.145,2); % ylocation

% Function for the attack
subplot('Position',[xxs(1) yys(1) wid hei]); 

testA=[zeros(21,1); 1; zeros(21,1)]'; % the days we want to plot
if(XU(3)+XU(6)*A~=0) % See if the function is present
    Y = ImpactAttack(testA,DB,DA,[0 0 0 0 0]); % evalautae function the function for the days (do not need lag here as we are interested in the function)
else
    Y=0.*testA; % the function is absent so set to zero
end
% Plot a bar graph since it is discrete
bar([-21:21],Y,'k','LineWidth',2);
% label for x-axis
xlabel('Day relative to the attack');
% lable for y-axis
ylabel({'Impact of water related',' attack on incidence'});
%remove box around figure
box off;
% adjust characteristics of axis
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[-21:2:21],'Ytick',[0:0.2:1],'Xminortick','on');
% x-axis limits
xlim([-21 21]);
% y-axis limits
ylim([0 1]);

%Function for conflict
subplot('Position',[xxs(2) yys(1) wid hei]);
Ct=linspace(0,450,501); % the values where the function will be evalauted
if(XU(4)+XU(6)*C~=0) % See if the function is present
    Y = ImpactConflict(Ct,K,n,CF); % evalautae function the function for the values
else
    Y=0.*Ct; % functino not present
end
% plot the function
plot(Ct,Y,'k','LineWidth',2);
% label for x-axis
xlabel('Amount of conflict');
% label for y-axis
ylabel({'Impact of general','conflict on incidence'});
%remove box around figure
box off;
% adjust characteristics of axis
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Ytick',[0:0.2:1]);
%y-axis limits
ylim([0 1]);

% plot rain fall function
subplot('Position',[xxs(1) yys(2) wid hei]);
Rt=linspace(0,10,101); % the values where the function will be evalauted
if(XU(5)+XU(6)*R~=0) % See if the function is present
    Y = ImpactRainfall(Rt,RF,rl,rh); % evalautae function the function for the values
else
    Y=0.*Rt; % functino not present
end
% plot the function
plot(Rt,Y,'k','LineWidth',2);
% label for x-axis
xlabel('Rainfall');
% label for y-axis
ylabel({'Impact of rainfall','on incidence'});
%remove box around figure
box off;
% adjust characteristics of axis
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Ytick',[0:0.2:1]);
%y-axis limits
ylim([0 1]);


%% Plot the national level incidence compared to the model predicted incidence

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

D=sum(IData,2); % Adds the incidence of all the governoantes together to provide the national level
M=sum(Yt,1); % sums across the different areas to produce the national level incidence
NW=length(D); % Determine the number of weeks that we have for the Governorates IData is NWxNG in size
scatter([1:NW],D,40,'k','filled'); hold on;
plot([1+max(tau):NW],M,'b','LineWidth',2); hold off;
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

%% Set up information for the projection
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,131); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,131);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
tA(tA>1)=1; % Find the weeks where there are more than one attack that occurs and assume that there is only a single attack
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

%% Run the projection
NW=length(NatIData); % number of week want to go out to
[Yt,Pt]= ModelProjection(beta,WI(GNZI,:),A,tA(GNZI,:),DB,DA,C,Ctv(GNZI,:),K,n,R,Rtv(GNZI,:),RF,rl,rh,tau,CF,ma,mc,mr,NW-length(WI(1,:))); % Run the projection
Mt=[Yt Pt]; % Combined the model fit with the model projection into one matrix

%% Produce the figure for the projection
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

%Plot national data in a scatter plot
scatter([1:NW],NatIData,40,'k','filled'); hold on;
%plot the model fit
M=sum(Yt,1); % Aggregate the incidence in the various areas for the fitted portion
NWF=length(WI(1,:)); % Number of weeks for the fitting
 plot([1+max(tau):NWF],M,'b','LineWidth',2);  % Plot the fitted portion of the mode
 
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

%% Plot data for the areas of interest
G=[2 3 8 11 21]; % Specify the number area of interest ranges from 1-22
GovFigureIncidence(G,Yt,WI(GNZI,:),GNZI)
