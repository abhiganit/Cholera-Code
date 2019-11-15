function [X] = CalcCovariates(WI,tA,DB,DA,Ctv,K,n,tau,maxtau,CF,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,Rtv,RF,r0,rm,beta)
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
% WaSH
X(1,:,:)=WPIN(:,(1+maxtau-tau(1)):(end-tau(1)));
X(2,:,:)=WPIN(:,(1+maxtau-tau(2)):(end-tau(2))).*ImpactAttack(tA,DB(1),DA(1),tau(2),maxtau).*WI(:,(1+maxtau-tau(2)):(end-tau(2)));%./(WI(:,(1+maxtau-tau(2)):(end-tau(2)))+KI(1));
X(3,:,:)=WPIN(:,(1+maxtau-tau(3)):(end-tau(3))).*(ImpactConflict(Ctv(:,(1+maxtau-tau(3)):(end-tau(3))),K(1),n(1),CF(1))).*WI(:,(1+maxtau-tau(3)):(end-tau(3)));%./(WI(:,(1+maxtau-tau(3)):(end-tau(3)))+KI(1));
X(4,:,:)=WPIN(:,(1+maxtau-tau(4)):(end-tau(4))).*ImpactConflict(Mt(:,(1+maxtau-tau(4)):(end-tau(4))),K(3),n(3),CF(1)).*WI(:,(1+maxtau-tau(4)):(end-tau(4)));%./(WI(:,(1+maxtau-tau(4)):(end-tau(4)))+KI(1));
X(5,:,:)=WPIN(:,(1+maxtau-tau(5)):(end-tau(5))).*max(Dieselt(:,(1+maxtau-tau(5)):(end-tau(5)))-KP(1),0).*WI(:,(1+maxtau-tau(5)):(end-tau(5)));%./(WI(:,(1+maxtau-tau(5)):(end-tau(5)))+KI(1));
%Food security
X(6,:,:)=FPIN(:,(1+maxtau-tau(6)):(end-tau(6)));
X(7,:,:)=FPIN(:,(1+maxtau-tau(7)):(end-tau(7))).*(ImpactConflict(Ctv(:,(1+maxtau-tau(7)):(end-tau(7))),K(2),n(2),CF(2))).*WI(:,(1+maxtau-tau(7)):(end-tau(7)));%./(WI(:,(1+maxtau-tau(7)):(end-tau(7)))+KI(2));
X(8,:,:)=FPIN(:,(1+maxtau-tau(8)):(end-tau(8))).*ImpactConflict(Mt(:,(1+maxtau-tau(8)):(end-tau(8))),K(4),n(4),CF(2)).*WI(:,(1+maxtau-tau(8)):(end-tau(8)));%./(WI(:,(1+maxtau-tau(8)):(end-tau(8)))+KI(2));
X(9,:,:)=FPIN(:,(1+maxtau-tau(9)):(end-tau(9))).*max(Dieselt(:,(1+maxtau-tau(9)):(end-tau(9)))-KP(1),0).*WI(:,(1+maxtau-tau(9)):(end-tau(9)));%./(WI(:,(1+maxtau-tau(9)):(end-tau(9)))+KI(2));
X(10,:,:)=FPIN(:,(1+maxtau-tau(10)):(end-tau(10))).*max(Wheatt(:,(1+maxtau-tau(10)):(end-tau(10)))-KP(2),0).*WI(:,(1+maxtau-tau(10)):(end-tau(10)));%./(WI(:,(1+maxtau-tau(10)):(end-tau(10)))+KI(2));

%Rainfall
if(RF(1)>=0)
    X(11,:,:)=WPIN(:,(1+maxtau-tau(11)):(end-tau(11))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(11)):(end-tau(11))),RF(1),r0).*WI(:,(1+maxtau-tau(11)):(end-tau(11)));%./(WI(:,(1+maxtau-tau(11)):(end-tau(11)))+KI(1));
end
if(RF(2)>=0)
    X(12,:,:)=WPIN(:,(1+maxtau-tau(12)):(end-tau(12))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(12)):(end-tau(12))),RF(2),rm).*WI(:,(1+maxtau-tau(12)):(end-tau(12)));
end
end

