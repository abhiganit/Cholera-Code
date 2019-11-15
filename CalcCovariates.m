function [X] = CalcCovariates(WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,a)
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

II=WI./(WI+KP); % Saturated incidence per captia

% Tau is assumed to be fixed here
PDG=P(:,(1+maxtau-tau(2)):(end-tau(2))).*WI(:,(1+maxtau-tau(2)):(end-tau(2))); % Populatino density
HFG=H(:,(1+maxtau-tau(3)):(end-tau(3))).*WI(:,(1+maxtau-tau(3)):(end-tau(3))); % Number of health sits per 10,000
It=(a.*WPIN(:,(1+maxtau-tau(4)):(end-tau(4)))+(1-a).*FPIN(:,(1+maxtau-tau(4)):(end-tau(4)))).*WI(:,(1+maxtau-tau(4)):(end-tau(4))); % Using the past incidence with a lag of tau weeks
Gt=(a.*WPIN(:,(1+maxtau-tau(5)):(end-tau(5)))+(1-a).*FPIN(:,(1+maxtau-tau(5)):(end-tau(5)))).*P(:,(1+maxtau-tau(5)):(end-tau(5))).*WI(:,(1+maxtau-tau(5)):(end-tau(5))); 
HWt=(a.*WPIN(:,(1+maxtau-tau(6)):(end-tau(6)))+(1-a).*FPIN(:,(1+maxtau-tau(6)):(end-tau(6)))).*H(:,(1+maxtau-tau(6)):(end-tau(6))).*WI(:,(1+maxtau-tau(6)):(end-tau(6)));
RCt=(a.*WPIN(:,(1+maxtau-tau(7)):(end-tau(7)))+(1-a).*FPIN(:,(1+maxtau-tau(7)):(end-tau(7)))).*repmat(RC,1,length(WI(1,(1+maxtau-tau(7)):(end-tau(7))))).*P(:,(1+maxtau-tau(7)):(end-tau(7))).*WI(:,(1+maxtau-tau(7)):(end-tau(7))); % Rebel control
% tau is estimated
ITAt=II(:,(1+maxtau-tau(8)):(end-tau(8))).*(WPIN(:,(1+maxtau-tau(8)):(end-tau(8)))).*ImpactAttack(tA,DB(1),DA(1),tau(8),maxtau); % Product of incidence and attacks 
ICt=II(:,(1+maxtau-tau(9)):(end-tau(9))).*(a.*WPIN(:,(1+maxtau-tau(9)):(end-tau(9)))+(1-a).*FPIN(:,(1+maxtau-tau(9)):(end-tau(9)))).*ImpactConflict(Ctv(:,(1+maxtau-tau(9)):(end-tau(9))),K(1),n(1),CF(1)); %Product of incidence and conflict
IAt=II(:,(1+maxtau-tau(10)):(end-tau(10))).*(a.*WPIN(:,(1+maxtau-tau(10)):(end-tau(10)))+(1-a).*FPIN(:,(1+maxtau-tau(10)):(end-tau(10)))).*ImpactAttack(Mt,DBE(1),DAE(1),tau(10),maxtau); % Product of incidence and attacks 
IRt=WPIN(:,(1+maxtau-tau(11)):(end-tau(11))).*II(:,(1+maxtau-tau(11)):(end-tau(11))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(11)):(end-tau(11))),RIF(1),rl(1)); %Product of incidence and rainfall
Rt=WPIN(:,(1+maxtau-tau(12)):(end-tau(12))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(12)):(end-tau(12))),RF(1),rh(1)).*WI(:,(1+maxtau-tau(12)):(end-tau(12))); %rainfall
CRIt=ImpactConflict(Ctv(:,(1+maxtau-tau(13)):(end-tau(13))),K(2),n(2),CF(2)).*ImpactRainfall(Rtv(:,(1+maxtau-tau(13)):(end-tau(13))),RIF(2),rl(2)).*II(:,(1+maxtau-tau(13)):(end-tau(13))); %Product of incidence and conflict
TARIt=ImpactAttack(tA,DB(2),DA(2),tau(14),maxtau).*ImpactRainfall(Rtv(:,(1+maxtau-tau(14)):(end-tau(14))),RIF(3),rl(3)).*II(:,(1+maxtau-tau(14)):(end-tau(14))); % Product of incidence and attacks 
ARIt=ImpactAttack(Mt,DBE(2),DAE(2),tau(15),maxtau).*ImpactRainfall(Rtv(:,(1+maxtau-tau(15)):(end-tau(15))),RIF(4),rl(4)).*II(:,(1+maxtau-tau(15)):(end-tau(15))); % Product of incidence and attacks 


WheatIt=FPIN(:,(1+maxtau-tau(16)):(end-tau(16))).*II(:,(1+maxtau-tau(16)):(end-tau(16))).*Wheatt(:,(1+maxtau-tau(16)):(end-tau(16))); % We can use the impact rainf all function for the threshold of the price
DieselIt=FPIN(:,(1+maxtau-tau(17)):(end-tau(17))).*II(:,(1+maxtau-tau(17)):(end-tau(17))).*Dieselt(:,(1+maxtau-tau(17)):(end-tau(17))); % We can use the impact rainf all function for the threshold of the price
WasHDieselIt=WPIN(:,(1+maxtau-tau(18)):(end-tau(18))).*II(:,(1+maxtau-tau(18)):(end-tau(18))).*Dieselt(:,(1+maxtau-tau(18)):(end-tau(18))); % We can use the impact rainf all function for the threshold of the price

X=zeros(18,length(PDG(:,1)),length(PDG(1,:)));

% Constant
X(1,:,:)=(a.*WPIN(:,(1+maxtau-tau(1)):(end-tau(1)))+(1-a).*FPIN(:,(1+maxtau-tau(1)):(end-tau(1)))).*P(:,(1+maxtau-tau(1)):(end-tau(1)));
X(2,:,:)=(a.*WPIN(:,(1+maxtau-tau(2)):(end-tau(2)))+(1-a).*FPIN(:,(1+maxtau-tau(2)):(end-tau(2)))).*repmat(RC,1,length(WI(1,(1+maxtau-tau(2)):(end-tau(2))))).*P(:,(1+maxtau-tau(2)):(end-tau(2))).*WI(:,(1+maxtau-tau(2)):(end-tau(2))); % Rebel control
X(3,:,:)=HFG;
X(4,:,:)=It;
X(5,:,:)=Gt;
X(6,:,:)=HWt;
X(7,:,:)=RCt;
X(8,:,:)=ITAt;
X(9,:,:)=ICt;
X(10,:,:)=IAt;
X(11,:,:)=IRt;
X(12,:,:)=Rt;
X(13,:,:)=CRIt;
X(14,:,:)=TARIt;
X(15,:,:)=ARIt;
X(16,:,:)=WheatIt;
X(17,:,:)=DieselIt;
X(18,:,:)=WasHDieselIt;
end

