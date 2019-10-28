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
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,~,GNZI,maxtau] = LoadYemenData;

NGS=floor(length(GNZI)*PDS);
Itemp=sum(WI(GNZI,:),2);
Itemp(RC(GNZI)==0)=0;
GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
% Find the top max
for ii=1:ceil(NGS/4)
   f=find(Itemp==max(Itemp)); % Find the maximum
   Itemp(f)=0; % set maximum to zero for it is no longer selected
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2);
Itemp(RC(GNZI)==1)=0;
% Find the top max
for ii=(ceil(NGS/4)+1):ceil(NGS/2)
   f=find(Itemp==max(Itemp)); % Find the maximum
   Itemp(f)=0; % set maximum to zero for it is no longer selected
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
Itemp(RC(GNZI)==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
% Find the minimum contributors
for ii=(ceil(NGS/2)+1):ceil(3.*NGS/4)
   f=find(Itemp==min(Itemp)); % Select minimum
   Itemp(f)=max(Itemp); % Set to maximum for not selected again
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
Itemp(RC(GNZI)==1)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
% Find the minimum contributors
for ii=(ceil(3.*NGS/4)+1):NGS
   f=find(Itemp==min(Itemp)); % Select minimum
   Itemp(f)=max(Itemp); % Set to maximum for not selected again
   GTF(ii)=f; % Record index
end
GTF=sort(GTF)'; % Gov. to used in the fitting of the model. We sort to keep order consistent with GNZI
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
lbps=[-32.*ones(1,length(XU)) zeros(1,8) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,5)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,8) ones(1,7) log10([ones(1,8) 20 3 20 3 120 120 120 120 120])]; % specify the upperbound for the parameters 

lb=[-32.*ones(1,length(XU)) ones(1,8) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,5) ]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=[ 5.*ones(1,length(XU)) 4.*ones(1,8) 2.*ones(1,7) log10([ones(1,8) 20 3 20 3 120 120 120 120 120])]; % specify the upperbound for the parameters 

IntC=[1:15]+length(XU);

%% Run the fitting algorithm
pars([1:8]+length(XU))=ceil(4.*pars([1:8]+length(XU)));
pars([9:15]+length(XU))=ceil(3.*pars([9:15]+length(XU)))-1;
options = optimoptions('ga','MaxGenerations',50000,'MaxStallGenerations',500,'UseParallel',true,'FunctionTolerance',10^(-8),'InitialPopulationMatrix',pars); %
optionsps = optimoptions('patternsearch','UseParallel',true,'Cache','on','SearchFcn','searchlhs','FunctionTolerance',10^(-10),'UseCompleteSearch',true);

[par] =ga(@(x)OFuncProGA(x,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),[]),length(pars),[],[],[],[],lb,ub,[],IntC,options); 
par(XU==0)=-30; % for the recursive componetnt
par([1:8]+length(XU))=par([1:8]+length(XU))./4-0.01; % such that they do not push on the boundary
par([9:15]+length(XU))= (par([9:15]+length(XU))+1)./3-0.01; % such that they do not push on the boundary
[par,fvalfit] =patternsearch(@(x)OFuncProPS(x,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),[]),par,[],[],[],[],lbps,ubps,[],optionsps); 
fvalfit=10.^fvalfit;
par(XU==0)=-30; % for the recursive componetnt
% Calculate the cross-validation error by calculating error for all areas
% w/ non-zero incidence and then subtracting the value of fvalfit

% Load data for the cross-Validation at the district level
CVE=(10.^OFuncProPS(par,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI,1:NW),RC(GNZI),H(GNZI,1:NW),WPIN(GNZI,1:NW),Mt(GNZI,1:NW),[]))-fvalfit;
end