function [par,fvalfit,CVE]=ProFittingGA(XU,PDS,CF,RF,pars)
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

[lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,pars,CF,RF);
options = optimoptions('ga','FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','MaxFunEvals',10^5,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);
opts= optimset('MaxIter',10^5,'MaxFunEvals',10^5,'TolFun',10^(-12),'TolX',10^(-12),'Display','off'); 


[par] =ga(@(x)OFuncProGA(x,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF),length(pars),[],[],[],[],lb,ub,[],IntC,options); 

[par] = ExpandPar(par,XU,CF,RF,1);
par(XU==0)=-30; % for the recursive componetnt

[~,~,~,~,~,par] = BoundsFitting(XU,par,CF,RF);

[par] =patternsearch(@(x)OFuncProPS(x,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF),par,[],[],[],[],lbps,ubps,[],optionsps); 
[par,fvalfit]=fmincon(@(x)OFuncProPS(x,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF),par,[],[],[],[],lbps,ubps,[],opts);
% compute cross validation error
CVE=(OFuncProPS(par,CF,WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),XU,maxtau,RC(GNZI(GTCV)),WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),RF));


[par] = ExpandPar(par,XU,CF,RF,1);
par(XU==0)=-30; % for the recursive componetnt
end