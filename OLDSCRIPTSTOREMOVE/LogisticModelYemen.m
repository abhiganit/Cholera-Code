function [Yt,X]= LogisticModelYemen(beta,tA,DB,DA,Ctv,K,n,tau,maxtau,CF,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,V1,V2,KV,dV,Rtv,RF,r0,WI,Pop,CI,DAR,w)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks
%===============================
% Input
%===============================

% WI - The incidence for the different areas for the given week
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
% maxtau- the maximum lag allowed for all the models such that they return
% the same amount of data
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % P - population density for the different govenorates
 % RC - Rebel control
 % H - The number of health facilities for the different govenorates
 % WPIN - People in Need WASH
 % FPIN - People in need food security
 % Mt - Shellings
 % Wheatt- Price of wheat
 % Dieselt - Proce of diesel
 % mln - the threshold for the price of disease
 % a - fraction for people in need of wash vs food secirity
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form
% X- the value of the covariates
% WPIN and FPIN - People in Need of WaSh and Food Security
% PDG - population density for the govenrate
% HFG- Health facilities in the govnorates
% It - Incidence last week
% Gt - Inicedence in the other govnorates
% HWt - Health care and people in need
% RCt - Rebel control
% ITAt - Product of incidence (Indicator) and targeted attacks 
% ICt - Product of incidence (Indicator) and conflict 
% IAt - Product of incidence (Indicator) and air and shelling attacks 
% IRt - Rainfall and wash (Indicator)
% Rt- Rainfall incidence and WASH 
% CRIt - Conflict, rianfall incidence, and WASH (Indicator)
% TARIt -Targeted Attacks, rianfall, incidence and WASH (Indicator)
% ARIt -Air/Shelling Attacks, rianfall, incidence and WASH (Indicator)
% Wheatt - proce of wheat and food insecurity (Indicator)
% Dieselt - proce of diesel and food insecurity (Indicator)
% WaSHDieselt- price of diesel and WASH (Indicator)

%% Input for regression model

[X] = CalcCovariates(tA,DB,DA,Ctv,K,n,tau,maxtau,CF,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,Rtv,RF,r0,WI,Pop,CI,DAR,w);
dV1=ImpactAttack(V1-V2,0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2,0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);
%% Output of regression model: the predicted weekly incidence of the model
Yt=zeros(size(squeeze(X(1,:,:))))';
for ii=1:length(beta)
    Yt=Yt+(1-EOVC).*beta(ii).*(squeeze(X(ii,:,:))');
end

end

