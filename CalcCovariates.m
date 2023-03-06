function [X] = CalcCovariates(tA,DB,DA,Ctv,K,n,tau,maxtau,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,Rtv,Temptv,r0,temp_0,WI,Pop,CI,DAR,w)
%CALCCOVARIATES Claculates the covariates for the regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% rh - THreshold for rainfall for the covariat of rainfall
% tau - the lag to use for the incidence and enviromental factors  
% maxtau- the maximum lag allowed for all the models such that they return
% the same amount of data
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


X=zeros(length(tau),length(WPIN(:,1)),length(WPIN(1,(1+maxtau-tau(1)):(end-tau(1)))));

%Temporal threshold for conflict
Ctv(:,1:n(1))=0;
Mt(:,1:n(2))=0;


% Targeted Attacks
for ii=1:maxtau
    X(ii,:,:)=WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii))).*(ImpactAttack(tA,DB(1),DA(1),tau(ii),maxtau)).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end
% General conflict
for ii=(maxtau+1):2*maxtau
    X(ii,:,:)=(w.*WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))+(1-w).*FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))).*(ImpactConflict(Ctv(:,(1+maxtau-tau(ii)):(end-tau(ii))),K(1))).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end
%Shelling and air attacks
for ii=(2.*maxtau+1):3*maxtau
    X(ii,:,:)=(w.*WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))+(1-w).*FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))).*(ImpactConflict(Mt(:,(1+maxtau-tau(ii)):(end-tau(ii))),K(2))).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end
% Diesel prices
for ii=(3.*maxtau+1):4*maxtau
    X(ii,:,:)=(w.*WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))+(1-w).*FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))).*max(Dieselt(:,(1+maxtau-tau(ii)):(end-tau(ii)))-KP(1),0).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end

% Wheat prices
for ii=(4.*maxtau+1):5*maxtau
    X(ii,:,:)=FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii))).*max(Wheatt(:,(1+maxtau-tau(ii)):(end-tau(ii)))-KP(2),0).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end

%Rainfall
for ii=(5.*maxtau+1):6*maxtau
    X(ii,:,:)=WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii))).*(ImpactRainfall(Rtv(:,(1+maxtau-tau(ii)):(end-tau(ii))),r0)).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end

%Temprature
for ii=(6.*maxtau+1):7*maxtau
    X(ii,:,:)=(w.*WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))+(1-w).*FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))).*(ImpactTemprature(Temptv(:,(1+maxtau-tau(ii)):(end-tau(ii))),temp_0)).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end

%Incidence per capita
for ii=(7.*maxtau+1):8*maxtau
    X(ii,:,:)=(w.*WPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))+(1-w).*FPIN(:,(1+maxtau-tau(ii)):(end-tau(ii)))).*WI(:,(1+maxtau-tau(ii)):(end-tau(ii))).*Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))./(Pop(:,(1+maxtau-tau(ii)):(end-tau(ii)))+DAR.*CI(:,(1+maxtau-tau(ii)):(end-tau(ii))));
end
end
