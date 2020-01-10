function [parfinal,fvalfit,CVE]=ProFittingKfold(XU,CF,RF,pars)
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
for mm=1:length(pars(:,1))
[lb,ub,lbps,ubps,IntC,parst(mm,:)] = BoundsFitting(XU,pars(mm,:),CF,maxtau);
end
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

C = combnk([1:21],3);

NN=length(C(:,1));
parfinal=zeros(NN,length(pars(1,:)));

optionsps = optimoptions('patternsearch','MaxFunEvals',10^5,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true,'Display','off');

options = optimoptions('ga','MaxGenerations',10^4,'FunctionTolerance',10^(-6),'InitialPopulationMatrix',parst,'Display','off','HybridFcn',{@patternsearch,optionsps});
CVE=zeros(NN,1);
fvalfit=zeros(NN,1);
%% Run the fitting algorithm
parfor kk=1:NN
    GTF=setdiff([1:21],C(kk,:));
    GTCV=C(kk,:);
    [par,fvalfit(kk)] =ga(@(x)OFuncProGA(x,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF,PopS(GNZI(GTF),1:NW),CI(GNZI(GTF),1:NW)),length(parst(1,:)),[],[],[],[],lb,ub,[],IntC,options); 
        
    % compute cross validation error
    CVE(kk)=(OFuncProPS(par,CF,WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),XU,maxtau,WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),RF,PopS(GNZI(GTCV),1:NW),CI(GNZI(GTCV),1:NW)));
    [parfinal(kk,:)] = ExpandPar(par,XU,CF,maxtau,1);
end

for kk=1:NN
    parfinal(kk,XU==0)=-30; % for the recursive componetnt
end
end