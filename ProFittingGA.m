function [par,fvalfit,RRS]=ProFittingGA(XU,PDS,pars)
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
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,GNZI,maxtau] = LoadYemenData;

NW=floor(153*PDS);
lbps=[-32.*ones(1,length(XU)) zeros(1,7) 0 0 0 -32.*ones(1,8)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,7) 1 1 1 log10([1 1 113 10 12 12 1 1])]; % specify the upperbound for the parameters 

lb=[-32.*ones(1,length(XU)) ones(1,7) 0 0 0 -32.*ones(1,8)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=[ 5.*ones(1,length(XU)) 4.*ones(1,7) 2 2 2 log10([1 1 113 10 12 12 1 1])]; % specify the upperbound for the parameters 

IntC=[1:10]+length(XU);

%% Run the fitting algorithm
pars([1:7]+length(XU))=ceil(4.*pars([1:7]+length(XU)));
pars([8:10]+length(XU))=ceil(3.*pars([8:10]+length(XU)))-1;
options = optimoptions('ga','MaxGenerations',10000,'MaxStallGenerations',100,'UseParallel',true,'FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','UseParallel',true,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);

[par] =ga(@(x)OFuncProGA(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI,1:NW),RC(GNZI),H(GNZI,1:NW),WPIN(GNZI,1:NW),Mt(GNZI,GNZI)),length(pars),[],[],[],[],lb,ub,[],IntC,options); 
par(XU==0)=-30; % for the recursive componetnt
par([1:7]+length(XU))=par([1:7]+length(XU))./4-0.01; % such that they do not push on the boundary
par([8:10]+length(XU))= (par([8:10]+length(XU))+1)./3-0.01; % such that they do not push on the boundary
[par,fvalfit] =patternsearch(@(x)OFuncProPS(x,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI,1:NW),RC(GNZI),H(GNZI,1:NW),WPIN(GNZI,1:NW),Mt(GNZI,GNZI)),par,[],[],[],[],lbps,ubps,[],optionsps); 
fvalfit=10.^fvalfit;
par(XU==0)=-30; % for the recursive componetnt
% Calculate the cross-validation error
RRS=10.^OFuncProPS(par,WI(GNZI,(NW+1-maxtau):end),tA(GNZI,(NW+1-maxtau):end),Ctv(GNZI,(NW+1-maxtau):end),Rtv(GNZI,(NW+1-maxtau):end),XU,maxtau,P(GNZI,(NW+1-maxtau):end),RC(GNZI),H(GNZI,(NW+1-maxtau):end),WPIN(GNZI,(NW+1-maxtau):end),Mt(GNZI,GNZI));
end