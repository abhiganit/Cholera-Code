function [par,fvalfit]=ProFittingGA(XU,CF,RF,pars)
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
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);



%% Run the fitting algorithm
for mm=1:length(pars(:,1))
[lb,ub,lbps,ubps,IntC,parst(mm,:)] = BoundsFitting(XU,pars(mm,:),CF,maxtau);
end
optionsps = optimoptions('patternsearch','MaxFunEvals',10^5,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-12),'UseCompleteSearch',true,'Display','off','PlotFcn', @psplotbestf);

options = optimoptions('ga','MaxGenerations',10^4,'FunctionTolerance',10^(-6),'InitialPopulationMatrix',parst,'Display','off','PlotFcn', @gaplotbestf,'HybridFcn',{@patternsearch,optionsps}); %
    [par,fvalfit] =ga(@(x)OFuncProGA(x,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW)),length(parst(1,:)),[],[],[],[],lb,ub,[],IntC,options); 
    
    [par] = ExpandPar(par,XU,CF,maxtau,1);
    
par(XU==0)=-30; % for the recursive componetnt
end