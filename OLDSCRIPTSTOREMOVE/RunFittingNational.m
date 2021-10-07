function RunFittingNational(XU,tau,AF,CF,RF,PF,PE,DT)
% Runs the fitting for the national level data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify what you want to use in in the regression model (1 = included, 0 =
%excluded) (1X9)
        % XU(1)- beta_0        
        % XU(2)- population density
        XU(2)=0; % Does not include population density as will not alter national (i.e. absorbed into beta_0)    
        % XU(3)- health facilities
        XU(3)=0; % Does not include health facilities as will not alter national trend (i.e. absorbed into beta_0)
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only
        % XU(9)- For Govnoerate level (external force of infection from
        % other govnoertaes)
        XU(9)=0; % Set to zero as this variable is not used in the national fitting
% tau -Specify a lag of all the factors that are integrated in the model
% (1X6)
    % tau(1) - Past incidence
    % tau(2) - Product of incidence and attacks
    % tau(3) - Product of incidence and conflict
    % tau(4) - Product of incidence and rainfall
    % tau(5) - Perciptiation only
    % tau(6) - Residual incidence (Not used here)
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
% PF - Plot only the fit if, otherwise not plot fit
% PE - Plot fucntions, otherwise not plot 
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
if isfile([pwd '\Tables\ModelFitSummaryNational.dat'])
     fmf = fopen([pwd '\Tables\ModelFitSummaryNational.dat'],'a'); % Writing to the file name
else
     fmf = fopen([pwd '\Tables\ModelFitSummaryNational.dat'],'w'); % Writing to the file name
     fprintf(fmf, 'Intercept \t Population \t HealthFacilities \t Incidence \t Attack \t Conflcit \t Rainfall \t Precipitation  \t AttackFunction \t ConflictFunction \t RainfallFunction \t ObjectiveValue \t AIC \n'); % The header for the table
end
%% Load the data
load('Yemen_National_Incidence.mat'); % Incidence data
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=NatIData'; % Transpose the data set such that the number of areas is the row
Ctv=sum(GLevelConflict(ProC,S,length(WI))); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
%% NEED TO ENSURE THAT THE ROWS ARE THE PROPER GOVERNORATES %%
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas
Rtv=sum(Rtv(:,1:length(NatIData))); % truncate the railfall data to the time period spceified
load('Attack_Yemen_Time_Location_Project_Forward.mat');
tA=sum(GLevelConflict(ProA,S,length(WI)));% % Send attack data (time, latitude, longatude) and shapefile of area wanted to catagorize
tA(tA>1)=1; % Find the weeks where there are more than one attack that occurs and assume that there is only a single attack

%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model
tau(tau>maxtau)=maxtau; % Set to ensure that the lag does not exceed the maximum lag specified
f=find(XU(2:end)==0); % Find the variables not included so the lag can be set to zero
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
problem = createOptimProblem('lsqnonlin','objective',@(x)OFunc(x,WI,tA,Ctv,Rtv,XU,tau,maxtau,AF,CF,RF,0,0),'x0',log10(rand(15,1)),'lb',lb,'ub',ub,'Aineq',[],'bineq',[],'Aeq',[],'beq',[],'options',opts);
ms = MultiStart('UseParallel','always'); % specifies that we run the algorithm in parallel
[par,fval] = run(ms,problem,NS); %starts at NS random initial points to thoroghly search the paramter space
% Evaluate the number of paramters that are being used in the estimation 
[k,beta,DB,DA,K,n,rl,rh]=RetParameter(par,XU,AF,CF,RF);

% Computing the statistical of the bets estimates for the regresison model
    f=find(XU==1); %  fiind components that are included
    be=beta(f); % increase index by one as we do not remove the inital beta_0
    [Yt,It,~,~,IAt,ICt,IRt,Rt,~]=LogisticModel(beta,WI,tA,DB,DA,Ctv,K,n,Rtv,RF,rl,rh,tau,maxtau,CF,0,0); % Returns the incidence in matrix form of size Ng X (NW-tau)
    resid=OFunc(par,WI,tA,Ctv,Rtv,XU,tau,maxtau,AF,CF,RF,0,0); % Returns the vector of residulas for each of the data points
    X=[]; % Construct the input vector based on the variables used in the fitting
    twotail=[]; % for the two-tail or single tail p-value (1 - uses two-tail)
    % If the constant is used
    if(XU(1)==1)
        X=[X ones(size(It(:)))]; % add the constant one to the input
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
    % If precipitation factors are used
    if(XU(8)==1)
        X=[X Rt(:)]; % add precipitation
        twotail=[twotail 1]; % Specify whether want two tail or not
    end
    
    % Calculate the statistics for for the coefficients of the model
    [SE,tStat,pValue] = LinRegressionStat(be,resid,k,X,twotail);   % We have specified the two-tail p-value as seen above
    
    % Constrcuts a table that is saved an can be viewed later
    f1 = fopen([pwd '\Tables\Table-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(2:end)>0)) '-National.dat'],'w'); % Writing to the file name
    CN=struct('N',{'Intercept ','Incidence ','Attack ','Conflcit ','Rainfall ','Precipitation '}); % The name of the variables that we can call
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
        readtable([pwd '\Tables\Table-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(2:end)>0)) '-National.dat']) % reads and prints the table to the command window
    end
    save([pwd '\Tables\Model-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(2:end)>0)) '-National.mat'],'k','beta','DB','DA','K','n','rl','rh','fval','tau','AIC'); % Saves the information for the specified model for can load and run model after
    
    if(DT~=0)
        fprintf('The objective value for the model fit: %5.4e \n',fval); % Outputs the AIC score
        fprintf('The AIC score for the model fit: %4.1f \n',AIC); % Outputs the AIC score
    end
fprintf(fmf, [ '%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t ' AFF((AF+1).*XU(3)+1).N '\t ' CFF((CF+1).*XU(4)+1).N ' \t ' RFF((RF+1).*XU(5)+1).N ' \t %8.0f \t %5.2f \n'],[XU(1:3) XU(4:end).*tau fval AIC]); % The header for the table
fclose(fmf);
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
    if(XU(3)~=0) % See if the function is present
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
    if(XU(4)~=0) % See if the function is present
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
    if(XU(5)~=0) % See if the function is present
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
    print(gcf,[pwd '\Figures\Functions-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(2:end)>0)) '-National.png'],'-dpng','-r600');
end

%% Plot the national level incidence compared to the model predicted incidence


if(PF~=0) % Plot fit 

    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

    D=WI; % Adds the incidence of all the governoantes together to provide the national level
    NW=length(D); % Determine the number of weeks that we have for the Governorates IData is NWxNG in size
    scatter([1:NW],D,40,'k','filled'); hold on;
    plot([1+maxtau:NW],Yt,'b','LineWidth',2); hold off;
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
    
    print(gcf,[pwd '\Figures\Fit-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF) '-AF=' num2str(AF) '-RF=' num2str(RF) '-tau=' num2str(tau(XU(2:end)>0)) '-National.png'],'-dpng','-r600');
end

end