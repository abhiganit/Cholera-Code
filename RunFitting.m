function RunFitting(XU,tau,AF,CF,RF,G,PF,PE,PP,DT)
% Runs the fitting for the specified criteria and saves files to folders
% for what is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify what you want to use in in the regression model (1 = included, 0 =
%excluded) (1X9)
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates
% tau -Specify a lag of all the factors that are integrated in the model
% (1X6)
    % tau(1) - Past incidence
    % tau(2) - Product of incidence and attacks
    % tau(3) - Product of incidence and conflict
    % tau(4) - Product of incidence and rainfall
    % tau(5) - Perciptiation only
    % tau(6) - Residual incidence
% AF -Specify the attack function to be used
        % AF=0 attack only has effect before;
        %AF=1 Attack has effect only after; 
        %AF=2; Attack has effect before and after
% CF - Specify the conflict function to be used
    % CF=0 linear effect; 
    % CF=1 Hill function with n=1;
    % CF=2; Full hill function; 
% RF- Specify the rainfall function to be used
    % RF=0 Increased incidence for low-rainfall; 
    % RF=1 increased incidence for high rainfall;
    % RF=2 increased incidence for high and low rain fall
% G- Specify the number area of interest ranges from 1-22
% PF - Plot only the fit if , otherwise not plot fit
% PE - Plot fucntions , otherwise not plot 
% PP - Plot fit and projection, otherwise not plot 
% DT - Display table/fitting projection results ro not 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves files to specified folders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Run alorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close all figures and clear all other data and clear command window
close all;
clc;
RFF=struct('N',{'None','Low','High','LowHigh'});
CFF=struct('N',{'None','Linear','Hill(n=1)','Hill'});
AFF=struct('N',{'None','Before','After','BeforeAfter'});
if isfile([pwd '\Tables\ModelFitSummary.dat'])
     fmf = fopen([pwd '\Tables\ModelFitSummary.dat'],'a'); % Writing to the file name
else
     fmf = fopen([pwd '\Tables\ModelFitSummary.dat'],'w'); % Writing to the file name
     fprintf(fmf, 'Intercept \t Incidence \t Attack \t Conflcit \t Rainfall \t Precipitation \t ExternalIncidence \t AttackFunction \t ConflictFunction \t RainfallFunction \t ObjectiveValue \t AIC \t MeanErrorProejection \t R2Projection \n'); % The header for the table
end
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

% Load population density
load('PopulationDensity_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
P=P';
% Load Health facility density
H=zeros(size(P));

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model
tau(tau>maxtau)=maxtau; % Set to ensure that the lag does not exceed the maximum lag specified
f=find(XU(4:end)==0); % Find the variables not included so the lag can be set to zero
tau(f)=0; % set the lag for components not included to zero



% Specify the lower bounds for the estimated parameters
% beta(1:9)=XU.*[x(1) 10.^(x(2:9))]
% 
% %Attack associated paramters
% DB=10.^x(10);
% DA=10.^x(11);
% % Conflict associated paramters
% K=10.^x(12);
% n=10.^x(13);
% % Rainfall assocaited paramters
% rl=10.^x(14);
% rh=10.^x(15);


lb=[-inf.*ones(1,15)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=log10([inf inf inf inf inf inf inf inf inf 1 1 inf inf 10 10]); % specify the upperbound for the parameters 


%% Run the fitting algorithm
NS=1000; % The number of initial starting points used in the fitting process
opts= optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-12),'TolX',10^(-12),'UseParallel',true,'Display','off');  % Specifies the conditions that will be used in the fitting processs
problem = createOptimProblem('lsqnonlin','objective',@(x)OFunc(x,WI(GNZI,:),tA(GNZI,:),Ctv(GNZI,:),Rtv(GNZI,:),XU,tau,maxtau,AF,CF,RF,P(GNZI),H(GNZI)),'x0',log10(rand(15,1)),'lb',lb,'ub',ub,'Aineq',[],'bineq',[],'Aeq',[],'beq',[],'options',opts);
ms = MultiStart('UseParallel','always'); % specifies that we run the algorithm in parallel
[par,fval] = run(ms,problem,NS); %starts at NS random initial points to thoroghly search the paramter space
% Evaluate the number of paramters that are being used in the estimation 
[k,beta,DB,DA,K,n,rl,rh]=RetParameter(par,XU,AF,CF,RF);

% Computing the statistical of the bets estimates for the regresison model
    f=find(XU==1); %  fiind components that are included
    be=beta(f); % increase index by one as we do not remove the inital beta_0
    [Yt,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt]=LogisticModel(beta,WI(GNZI,:),tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,Rtv(GNZI,:),RF,rl,rh,tau,maxtau,CF,P(GNZI),H(GNZI)); % Returns the incidence in matrix form of size Ng X (NW-tau)
    resid=OFunc(par,WI(GNZI,:),tA(GNZI,:),Ctv(GNZI,:),Rtv(GNZI,:),XU,tau,maxtau,AF,CF,RF,P(GNZI),H(GNZI)); % Returns the vector of residulas for each of the data points
    X=[]; % Construct the input vector based on the variables used in the fitting
    twotail=[]; % for the two-tail or single tail p-value (1 - uses two-tail)
    % If the constant is used
    if(XU(1)==1)
        X=[X ones(size(It(:)))]; % add the constant one to the input
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    if(XU(2)==1) % Population density
        X=[X PDG(:)]; % add the constant one to the input
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    if(XU(3)==1)  % Health facilities 
        X=[X HFG(:)]; % add the constant one to the input
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    %If past incidence is used
    if(XU(4)==1)
        X=[X It(:)]; % add the past incidence
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    %If attacks are used
    if(XU(5)==1)
        X=[X IAt(:) ]; % add product of past incidence and attacks
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If confclit is used
    if(XU(6)==1)
        X=[X ICt(:)]; % add product of past incidence and conflict
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If rainfall is used
    if(XU(7)==1)
        X=[X IRt(:)]; % add product of past incidence and rainfall
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % If percipitation is used
    if(XU(8)==1)
        X=[X Rt(:)]; % add percipitation is used
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    
    % If incidence of other govneroates are used
    if(XU(9)==1)
        X=[X Gt(:)]; % add product of past incidence and environmental factors
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    % Calculate the statistics for for the coefficients of the model
    [SE,tStat,pValue] = LinRegressionStat(be,resid,k,X,twotail);   % We have specified the two-tail p-value as seen above
    
    % Constrcuts a table that is saved an can be viewed later
    f1 = fopen([pwd '\Tables\Table-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.dat'],'w'); % Writing to the file name
    CN=struct('N',{'Intercept ','Population','Health','Incidence ','Attack ','Conflcit ','Rainfall ','Precipitation ','External '}); % The name of the variables that we can call
    fprintf(f1, 'Coefficient \t Value \t Stand._Error \t t-Stat \t p-value \n'); % The header for the table
    for ii=1:length(XU) % Go through all variables to see which are included
        if(XU(ii)~=0) % print on the variable that are included in the model
            fprintf(f1, [CN(ii).N ' \t %4.3f \t %4.3f \t %4.3f \t %3.2E \n'],[be(sum(XU(1:ii))) SE(sum(XU(1:ii))) tStat(sum(XU(1:ii))) pValue(sum(XU(1:ii)))]); % writes the information to the file
        end
    end
    fclose(f1); % closes the file
    % Calculate the AIC score for the model fit
    nData=length(Yt(:)); % The number of data points available for the model fitting
    AIC = AICScore(k,nData,fval); % Calculate the AIC score
    if(DT~=0)
        readtable([pwd '\Tables\Table-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.dat']) % reads and prints the table to the command window
    end
    save([pwd '\Tables\Model-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.mat'],'k','beta','DB','DA','K','n','rl','rh','fval','tau','AIC'); % Saves the information for the specified model for can load and run model after
    
    if(DT~=0)
        fprintf('The objective value for the model fit: %5.4e \n',fval); % Outputs the AIC score
        fprintf('The AIC score for the model fit: %4.1f \n',AIC); % Outputs the AIC score
    end
%% Plot the enviromental functions for the model

if(PE~=0)
    figure('units','normalized','outerposition',[0 0 1 1]);

    % Dimensions for the subpanels 
    hei=0.23; % heigh of figure
    wid=0.84; % width of figure
    xxs=0.075; % x-location
    yys=linspace(0.97-hei,0.08,3); % ylocation

    % Function for the attack
    subplot('Position',[xxs(1) yys(1) wid hei]); 

    testA=[zeros(21,1); 1; zeros(21,1)]'; % the days we want to plot
    if(XU(5)~=0) % See if the function is present
        Y = ImpactAttack(testA,DB,DA,[0 0 0 0 0],0); % evalautae function the function for the days (do not need lag here as we are interested in the function)
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
    subplot('Position',[xxs(1) yys(2) wid hei]);
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
    set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'Ytick',[0:0.2:1]);
    %y-axis limits
    ylim([0 1]);

    % plot rain fall function
    subplot('Position',[xxs(1) yys(3) wid hei]);
    Rt=linspace(0,10,101); % the values where the function will be evalauted
    if(XU(7)~=0) % See if the function is present
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
    print(gcf,[pwd '\Figures\Functions-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.png'],'-dpng','-r600');
end

%% Plot the national level incidence compared to the model predicted incidence


if(PF~=0) % Plot fit 

    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

    D=sum(IData,2); % Adds the incidence of all the governoantes together to provide the national level
    M=sum(Yt,1); % sums across the different areas to produce the national level incidence
    NW=length(D); % Determine the number of weeks that we have for the Governorates IData is NWxNG in size
    scatter([1:NW],D,40,'k','filled'); hold on;
    plot([1+maxtau:NW],M,'b','LineWidth',2); hold off;
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
    
    print(gcf,[pwd '\Figures\Fit-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.png'],'-dpng','-r600');
end

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
[Yt,Pt]= ModelProjection(beta,WI(GNZI,:),tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,Rtv(GNZI,:),RF,rl,rh,tau,maxtau,CF,P(GNZI),H(GNZI),NW-length(WI(1,:))); % Run the projection
Mt=[Yt Pt]; % Combined the model fit with the model projection into one matrix
%% Compute the R^2 value of the model projection
RPro=corr(sum(Pt,1)',NatIData((length(WI(1,:))+1):end)).^2; % The R^2 value for the projection component of the model fit
RFit=corr(sum(Yt,1)',NatIData(1+maxtau:length(WI(1,:)))).^2; % The R^2 value for the projection component of the model fit

if(DT~=0)
    fprintf('R^2 value for the model fit: %3.2f \n',RFit);
    fprintf('R^2 value for the model projection: %3.2f \n',RPro);
end
MSEP=mean((sum(Pt,1)'-NatIData((length(WI(1,:))+1):end)).^2);
if(DT~=0)
    fprintf('Mean Error in model projection: %3.0f \n',MSEP);
end

save([pwd '\Tables\Model-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.mat'],'RPro','RFit','MSEP','-append'); % Saves the information for the specified model for can load and run model after

% Save to model fit summary

fprintf(fmf, [ '%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t ' AFF((AF+1).*XU(5)+1).N '\t ' CFF((CF+1).*XU(6)+1).N ' \t ' RFF((RF+1).*XU(7)+1).N ' \t %8.0f \t %5.2f \t %8.0f \t %3.2f \n'],[XU(1:3) XU(4:end).*tau fval AIC MSEP RPro]); % The header for the table
fclose(fmf);
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
    print(gcf,[pwd '\Figures\FitandProjection-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(4:end)>0)) '.png'],'-dpng','-r600');
end

%% Plot data for the areas of interest
    if(~isempty(G))
         GovFigureIncidence(G,Yt,WI(GNZI,:),GNZI)
    end
end