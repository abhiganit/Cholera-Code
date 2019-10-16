function [par,fvalfit,CVE]=ProFittingGA(XU,PDS,Gov,pars)
% Runs the fitting for the specified criteria and saves files to folders
% for what is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify the model being fit    
        % XU(1)- beta_0
        % XU(2) - population density for the govenrate
        % XU(3)- Health facilities in the govnorates
        % XU(4) - Incidence last week
        % XU(5) - Inicedence in the other govnorates
        % XU(6) - Internally displaced people
        % XU(7) - Rebel control
        % XU(8) - Product of incidence and attacks 
        % XU(9) - Product of incidence and conflict 
        % XU(10) - Product of cumulative attacks incidence and rainfall 
        % XU(11)- cumulative attacks and Rainfall
        % XU(12) - Conflict, rianfall incidence
        % XU(13) - Attack, rianfall, incidence
% PDS - Percentage of governorates to be used in the fitting of the model (0<=PDS<=1)
% pars - starting point for the pattern search algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves files to specified folders
% par - the vector of log_10 estimated parameters from the fitting
% fval - the residual sum of squares
% CVE - the cross-validation error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Run alorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close all figures and clear all other data and clear command window
close all;
clc;

%% Load the data
[WI,Ctv,tA,Rtv,Mt,maxtau] = LoadYemenData;
NW=floor(153*PDS);
lbps=[-32.*ones(1,length(XU)) zeros(1,8) zeros(1,6) -32.*ones(1,6) -32.*ones(1,4) -32.*ones(1,4)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,8) ones(1,6) log10([ones(1,6) 20 3 20 3 120 120 120 120])]; % specify the upperbound for the parameters 

lb=[-32.*ones(1,length(XU)) ones(1,8) zeros(1,6) -32.*ones(1,6) -32.*ones(1,4) -32.*ones(1,4) ]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=[ 5.*ones(1,length(XU)) 4.*ones(1,8) 2.*ones(1,6) log10([ones(1,6) 20 3 20 3 120 120 120 120])]; % specify the upperbound for the parameters 

IntC=[1:14]+length(XU);

%% Run the fitting algorithm
pars([1:8]+length(XU))=ceil(4.*pars([1:8]+length(XU)));
pars([9:14]+length(XU))=ceil(3.*pars([9:14]+length(XU)))-1;
options = optimoptions('ga','MaxGenerations',50000,'MaxStallGenerations',500,'UseParallel',true,'FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','UseParallel',true,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);

[par] =ga(@(x)OFuncProGA(x,WI(Gov,1:NW),tA(Gov,1:NW),Ctv(Gov,1:NW),Rtv(Gov,1:NW),XU,maxtau,Mt(Gov,1:NW)),length(pars),[],[],[],[],lb,ub,[],IntC,options); 
par(XU==0)=-30; % for the recursive componetnt
par([1:8]+length(XU))=par([1:8]+length(XU))./4-0.01; % such that they do not push on the boundary
par([9:14]+length(XU))= (par([9:14]+length(XU))+1)./3-0.01; % such that they do not push on the boundary
[par,fvalfit] =patternsearch(@(x)OFuncProPS(x,WI(Gov,1:NW),tA(Gov,1:NW),Ctv(Gov,1:NW),Rtv(Gov,1:NW),XU,maxtau,Mt(Gov,1:NW)),par,[],[],[],[],lbps,ubps,[],optionsps); 
par(XU==0)=-30; % for the recursive componetnt
% Calculate the cross-validation error by calculating error for all areas
% w/ non-zero incidence and then subtracting the value of fvalfit
CVE=(OFuncProPS(par,WI(Gov,1:end),tA(Gov,1:end),Ctv(Gov,1:end),Rtv(Gov,1:end),XU,maxtau,Mt(Gov,1:end)))-fvalfit;
end