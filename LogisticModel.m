function [Yt,X]= LogisticModel(beta,WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,Mt,IDPt)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks
%===============================
% Input
%===============================
% t - the single time of interest
% WI - The incidence for the different areas for the given week
% beta- the coefficients for the regression model
% tA- A  vector indicating the time of the attacks 
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% DBE - the number of days before the attack that environemnt is affected
% DAE - the number of days after the attack environment affected
% Ct - the conflict at time t
% K- the saturation function
% n - the hill coefficient
% Rt - the rainfall at time t
% RIF - the indication of what function will be used
    % RIF=0 the increase in incidence when rainfall is low
    % RIF=1 the increase in incidence when rainfall is high
    % RIF=2 the increase in incidence when rainfall is low and high
% rl - THreshold for rainfall for the covariat of rainfall*incidence
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
% rh - THreshold for rainfall for the covariat of rainfall
% tau - the lag to use for the incidence and enviromental factors    
    %tau(1) - population density incidence
    % tau(2) - health zone incidence
    % tau(3) - Past incidence
    % tau(4) - Product of incidence and attacks
    % tau(5) - Product of incidence and conflict
    % tau(6) - Product of incidence and rainfall
    % tau(7) - Perciptiation only
    % tau(8) - Residual incidence
    % tau(9) - Attacks only
% maxtau- the maximum lag allowed for all the models such that they return
% the same amount of data
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % P - population density for the different govenorates
 % RC - Rebel control
 % H - The number of health facilities for the different govenorates
 % WPIN - Wash people in need
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form
% X- the value of the covariates
% Constnat - ones 
% PDG - population density for the govenrate
% HFG- Health facilities in the govnorates
% It - Incidence last week
% IAt - Product of incidence and attacks 
% ICt - Product of incidence and conflict 
% IRt - Product of cumulative attacks incidence and rainfall 
% Rt- cumulative attacks and Rainfall
% Gt - Inicedence in the other govnorates
% At- Attack only
% RCt - Rebel control
% WPt - WASH Status
% WPIt - Wash status and incidence

%% Input for regression model

[X] = CalcCovariates(WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,Mt,IDPt);

%% Output of regression model: the predicted weekly incidence of the model
Yt=zeros(size(squeeze(X(1,:,:))));
for ii=1:length(beta)
    Yt=Yt+beta(ii).*squeeze(X(ii,:,:));
end

end

