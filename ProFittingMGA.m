function [par,fvalfit,RRS]=ProFittingMGA(XU,PDS,pars)
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
% RRS - Cross-validation error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Run alorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,GNZI,maxtau] = LoadYemenData;

NW=floor(153*PDS); % Determine the number of points to be used in the fitting

% Bounds for the fitting
lbps=[-32.*ones(1,length(XU)) 10^(-6).*ones(1,10) -32.*ones(1,8)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,7) 1 1 1 log10([1 1 113 10 12 12 1 1])]; % specify the upperbound for the parameters 


%% Run the fitting algorithm
options = optimoptions('gamultiobj','MaxGenerations',50000,'MaxStallGenerations',100,'UseParallel',true,'InitialPopulationMatrix',pars); %
[par] =gamultiobj(@(x)MOFuncProGA(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI),Mt(GNZI,GNZI)),length(pars),[],[],[],[],lbps,ubps,[],options); 
options  = optimoptions('paretosearch','UseParallel',true,'ParetoSetChangeTolerance',10^(-4),'InitialPoints',par); %
[par,fval] =paretosearch(@(x)MOFuncProGA(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI),Mt(GNZI,GNZI)),length(pars),[],[],[],[],lbps,ubps,[],options); 
fval=sum(10.^fval,2); % Transform fval 
    
%% Calculate the cross validation error
RRSv=zeros(length(fval),1);
for mm=1:length(RRSv) % the minus tau is because we have to use the data points in the past to predict the incidence for the cross validation
    RRSv(mm)=10.^OFuncProPS(par(mm,:),WI(GNZI,(NW+1-maxtau):end),tA(GNZI,(NW+1-maxtau):end),Ctv(GNZI,(NW+1-maxtau):end),Rtv(GNZI,(NW+1-maxtau):end),XU,maxtau,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI),Mt(GNZI,GNZI));
end
gg=find(fval==min(fval),1); % find the pareto optimal soulation with the lowest error as we use this in the f-test
RRS=RRSv(gg); % record cross validation error
fvalfit=fval(gg); % record lowest error
par=par(gg,:); % record paramter values
par(XU==0)=-30; % for the recursive componetnt to ensure that they are small when fitting later models 

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
save(['ModelSelection-MultiObj-PercentDataSet=' num2str(PDS*100) '-XU=' num2str(XU*((2.^[0:(length(XU)-1)])')) '-CF=' num2str(CF*XU(6)) '-RIF=' num2str(RIF*XU(7)) '-PF=' num2str(RF*XU(8)) '.mat'],'k','DBE','DAE','beta','DB','DA','K','n','rl','rh','fvalfit','tau','par','RRS'); % Saves the information for the specified model for can load and run model after
end