function [par,fvalfit,RRS]=ProFittingGA(XU,PDS,G,PE,PP,DT,pars)
% Runs the fitting for the specified criteria and saves files to folders
% for what is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify what you want to use in in the regression model (1 = included, 0 =
%excluded) (1X11)
        % Specify XU the model being fit
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates
        % XU(10) - Attacks only
        % XU(11) - Rebel control
% PDS - Percentage of the data set to be used in the fitting of the model (0<=PDS<=1)
% G- Specify the number area of interest ranges from 1-22
% PF - Plot only the fit if , otherwise not plot fit
% PE - Plot fucntions , otherwise not plot 
% PP - Plot fit and projection, otherwise not plot 
% DT - Display table/fitting projection results ro not 
% pars - starting point for the pattern search algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves files to specified folders
% par - the vector of log_10 estimated parameters from the fitting
% fval - the residual sum of squares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Run alorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close all figures and clear all other data and clear command window
close all;
clc;

%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
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
% WASH People in Need Desnity
load('WASH_PeopleinNeed_Density_Yemen.mat');
WPIN=log(1+WPIN)';
%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model



% Specify the lower bounds for the estimated parameters
% beta=[10.^x(1:11)'].*XU;
% 
% %Attack associated paramters
% DB=10.^x(12);
% DA=10.^x(13);
% % Conflict associated paramters
% K=10.^x(14);
% n=10.^x(15);
% % Rainfall assocaited paramters
% rl=10.^x(16);
% rh=10.^x(17);
% DAE=10.^x(18);
NW=floor(153*PDS);
lbps=[-16.*ones(1,11) zeros(1,5) 0 0 0 -16.*ones(1,8)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,11) ones(1,5) 1 1 1 log10([1 1 113 10 12 12 1 1])]; % specify the upperbound for the parameters 

lb=[-16.*ones(1,11) ones(1,5) 0 0 0 -16.*ones(1,8)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=[ 5.*ones(1,11) 4.*ones(1,5) 2 2 2 log10([1 1 113 10 12 12 1 1])]; % specify the upperbound for the parameters 

IntC=[12:19];

%% Run the fitting algorithm
pars(12:16)=ceil(4.*pars(12:16));
pars(17:19)=ceil(3.*pars(17:19))-1;
options = optimoptions('ga','MaxGenerations',10000,'MaxStallGenerations',100,'UseParallel',true,'FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','UseParallel',true,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);

[par] =ga(@(x)OFuncProGA(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI)),length(pars),[],[],[],[],lb,ub,[],IntC,options); 
par(XU==0)=-15.99; % for the recursive componetnt
par(12:16)=par(12:16)./4-0.01; % such that they do not push on the boundary
par(17:19)= (par(17:19)+1)./3-0.01; % such that they do not push on the boundary
[par,fvalfit] =patternsearch(@(x)OFuncProPS(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI)),par,[],[],[],[],lbps,ubps,[],optionsps); 
fvalfit=10.^fvalfit;
par(XU==0)=-15.99; % for the recursive componetnt
% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);

    % Calculate the AIC score for the model fit
    save([pwd '\Tables\BackSelectModel-PercentDataSet=' num2str(PDS*100) '-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF*XU(6)) '-RIF=' num2str(RIF*XU(7)) '-PF=' num2str(RF*XU(8)) '-tau=' num2str(tau(XU(2:end)>0)) '.mat'],'k','beta','DB','DA','K','n','rl','rh','fvalfit','tau','par'); % Saves the information for the specified model for can load and run model after
    
%% Plot the enviromental functions for the model

if(PE~=0)
    figure('units','normalized','outerposition',[0 0 1 1]);

    % Dimensions for the subpanels 
    hei=0.39; % heigh of figure
    wid=0.42; % width of figure
    xxs=linspace(0.075,0.97-wid,2); % x-location
    yys=linspace(0.97-hei,0.08,2); % ylocation

    % Function for the attack
    subplot('Position',[xxs(1) yys(1) wid hei]); 

    testA=[zeros(21,1); 1; zeros(21,1)]'; % the days we want to plot
    if(XU(5)~=0) % See if the function is present
        Y1 = ImpactAttack(testA,DB,DA,0,0); % evalautae function the function for the days (do not need lag here as we are interested in the function)
    else
        Y1=0.*testA; % the function is absent so set to zero
    end
    if(XU(10)~=0) % See if the function is present
        Y2 = ImpactAttack(testA,0,DAE,0,0); % evalautae function the function for the days (do not need lag here as we are interested in the function)
    else
        Y2=0.*testA; % the function is absent so set to zero
    end
    Y=[Y1;Y2]';
    % Plot a bar graph since it is discrete
    bar([-21:21],Y,'LineWidth',2);
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
    if(XU(6)~=0) % See if the function is present
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
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Ytick',[0:25:450]);
    %y-axis limits
    ylim([0 450]);

    % plot rain fall function
    subplot('Position',[xxs(1) yys(2) wid hei]);
    Rt=linspace(0,16,101); % the values where the function will be evalauted
    if(XU(7)~=0) % See if the function is present
        Y = ImpactRainfall(Rt,RIF,rl); % evalautae function the function for the values
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
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
    %y-axis limits
    ylim([0 max(Y)+1]);
    xlim([0 16]);
    
    % plot rain fall function
    subplot('Position',[xxs(2) yys(2) wid hei]);
    Rt=linspace(0,16,101); % the values where the function will be evalauted
    if(XU(8)~=0) % See if the function is present
        Y = ImpactRainfall(Rt,RF,rh); % evalautae function the function for the values
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
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',16);
    %y-axis limits
    ylim([0 max(Y)+1]);    
    xlim([0 16]);
    print(gcf,[pwd '\Figures\Functions-PercentDataSet=' num2str(PDS*100) '-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF*XU(6)) '-RIF=' num2str(RIF*XU(7)) '-PF=' num2str(RF*XU(8)) '-tau=' num2str(tau(XU(2:end)>0)) '.png'],'-dpng','-r600');
end

%% Set up information for the projection

RRS=10.^OFuncProPS(par,WI(GNZI,:),tA(GNZI,:),Ctv(GNZI,:),Rtv(GNZI,:),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI));

save([pwd '\Tables\BackSelectModel-PercentDataSet=' num2str(PDS*100) '-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF*XU(6)) '-RIF=' num2str(RIF*XU(7)) '-PF=' num2str(RF*XU(8)) '-tau=' num2str(tau(XU(2:end)>0)) '.mat'],'RRS','-append'); % Saves the information for the specified model for can load and run model after


%% Produce the figure for the projection

if(PP~=0)
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
    print(gcf,[pwd '\Figures\ProjectionFit-PercentDataSet=' num2str(PDS*100) '-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF*XU(6)) '-RIF=' num2str(RIF*XU(7)) '-PF=' num2str(RF*XU(8)) '-tau=' num2str(tau(XU(2:end)>0)) '.png'],'-dpng','-r600');
end

%% Plot data for the areas of interest
    if(~isempty(G))
         GovFigureIncidence(G,Yt,WI(GNZI,:),GNZI)
    end
end