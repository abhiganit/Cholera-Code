function [par,fvalfit,CVE]=ProFittingGA(XU,PDS,pars)
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
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,PDS);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);



%% Run the fitting algorithm
pars([1:(length(XU)-7)]+length(XU))=ceil(4.*pars([1:(length(XU)-7)]+length(XU)));
pars([1:7]+2.*length(XU)-7)=ceil(3.*pars([1:7]+2.*length(XU)-7))-1;

[lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,pars);
options = optimoptions('ga','MaxGenerations',25000,'MaxStallGenerations',500,'UseParallel',false,'FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','UseParallel',false,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);

[par] =ga(@(x)OFuncProGA(x,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW)),length(pars),[],[],[],[],lb,ub,[],IntC,options); 

[par] = ExpandPar(par,XU,1);
par(XU==0)=-30; % for the recursive componetnt
par([1:(length(XU)-7)]+length(XU))=par([1:(length(XU)-7)]+length(XU))./4-0.01; % such that they do not push on the boundary
par([1:7]+2.*length(XU)-7)= (par([1:7]+2.*length(XU)-7)+1)./3-0.01; % such that they do not push on the boundary

[~,~,~,~,~,par] = BoundsFitting(XU,par);

[par,fvalfit] =patternsearch(@(x)OFuncProPS(x,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW)),par,[],[],[],[],lbps,ubps,[],optionsps); 

% Calculate the cross-validation error by calculating error for all areas
% w/ non-zero incidence and then subtracting the value of fvalfit

% Load data for the cross-Validation at the district level
CVE=(OFuncProPS(par,WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),XU,maxtau,P(GNZI(GTCV),1:NW),RC(GNZI(GTCV)),H(GNZI(GTCV),1:NW),WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW)));


[par] = ExpandPar(par,XU,1);
par(XU==0)=-30; % for the recursive componetnt
end