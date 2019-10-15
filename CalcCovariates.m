function [X] = CalcCovariates(WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,Mt,IDPt)
%CALCCOVARIATES Claculates the covariates for the regression model
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
    %tau(2) - population density incidence
    % tau(3) - health zone incidence
    % tau(4) - Past incidence
    % tau(5) - Incidence other governorates
    % tau(6) - Internally displaced
    % tau(7) - Rebel control
    % tau(8) - Product of incidence and attacks 
    % tau(9) - Product of incidence and conflict 
    % tau(10) - Product of cumulative attacks incidence and rainfall 
    % tau(11)- cumulative attacks and Rainfall
    % tau(12) - Conflict, rianfall incidence
    % tau(13) - Attack, rianfall, incidence
% maxtau- the maximum lag allowed for all the models such that they return
% the same amount of data
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % P - population density for the different govenorates
 % RC - Rebel control
 % H - The number of health facilities for the different govenorates
 % Mt - Contact matrix for the watercourse
%=================================
% Output
%=================================
% X- the value of the covariates
% Constnat - ones 
% PDG - population density for the govenrate
% HFG- Health facilities in the govnorates
% It - Incidence last week
% Gt - Inicedence in the other govnorates
% IDP - Internally displaced people
% RCt - Rebel control
% IAt - Product of incidence and attacks 
% ICt - Product of incidence and conflict 
% IRt - Product of cumulative attacks incidence and rainfall 
% Rt- cumulative attacks and Rainfall
% WPt - Conflict, rianfall incidence
% WPIt - Attack, rianfall, incidence

% Tau is assumed to be fixed here
PDG=P(:,(1+maxtau-tau(1)):(end-tau(1))).*WI(:,(1+maxtau-tau(1)):(end-tau(1))); % Populatino density
HFG=H(:,(1+maxtau-tau(2)):(end-tau(2))).*WI(:,(1+maxtau-tau(2)):(end-tau(2))); % Number of health sits per 10,000
It=WPIN(:,(1+maxtau-tau(3)):(end-tau(3))).*WI(:,(1+maxtau-tau(3)):(end-tau(3))); % Using the past incidence with a lag of tau weeks
Gt=WPIN(:,(1+maxtau-tau(4)):(end-tau(4))).*P(:,(1+maxtau-tau(4)):(end-tau(4))).*WI(:,(1+maxtau-tau(4)):(end-tau(4))); 
HWt=H(:,(1+maxtau-tau(5)):(end-tau(5))).*WPIN(:,(1+maxtau-tau(5)):(end-tau(5))).*WI(:,(1+maxtau-tau(5)):(end-tau(5)));
IDP=IDPt(:,(1+maxtau-tau(6)):(end-tau(6)));
RCt=repmat(RC,1,length(WI(1,(1+maxtau-tau(7)):(end-tau(7))))).*WI(:,(1+maxtau-tau(7)):(end-tau(7))); % Rebel control
% tau is estimated
ITAt=WI(:,(1+maxtau-tau(8)):(end-tau(8))).*ImpactAttack(tA,DB(1),DA(1),tau(8),maxtau); % Product of incidence and attacks 
ICt=WI(:,(1+maxtau-tau(9)):(end-tau(9))).*ImpactConflict(Ctv(:,(1+maxtau-tau(9)):(end-tau(9))),K(1),n(1),CF(1)); %Product of incidence and conflict
IAt=WI(:,(1+maxtau-tau(10)):(end-tau(10))).*ImpactAttack(Mt,DBE(1),DAE(1),tau(10),maxtau); % Product of incidence and attacks 
IRt=WPIN(:,(1+maxtau-tau(11)):(end-tau(11))).*WI(:,(1+maxtau-tau(11)):(end-tau(11))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(11)):(end-tau(11))),RIF(1),rl(1)); %Product of incidence and rainfall
Rt=WPIN(:,(1+maxtau-tau(12)):(end-tau(12))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(12)):(end-tau(12))),RF(1),rh(1)); %rainfall
CRIt=ImpactConflict(Ctv(:,(1+maxtau-tau(13)):(end-tau(13))),K(2),n(2),CF(2)).*ImpactRainfall(Rtv(:,(1+maxtau-tau(13)):(end-tau(13))),RIF(2),rl(2)).*WI(:,(1+maxtau-tau(13)):(end-tau(13))); %Product of incidence and conflict
TARIt=ImpactAttack(tA,DB(2),DA(2),tau(14),maxtau).*ImpactRainfall(Rtv(:,(1+maxtau-tau(14)):(end-tau(14))),RIF(3),rl(3)).*WI(:,(1+maxtau-tau(14)):(end-tau(14))); % Product of incidence and attacks 
ARIt=ImpactAttack(Mt,DBE(2),DAE(2),tau(15),maxtau).*ImpactRainfall(Rtv(:,(1+maxtau-tau(15)):(end-tau(15))),RIF(4),rl(4)).*WI(:,(1+maxtau-tau(15)):(end-tau(15))); % Product of incidence and attacks 

X=zeros(15,length(PDG(:,1)),length(PDG(1,:)));

% Constant
X(1,:,:)=PDG;
X(2,:,:)=HFG;
X(3,:,:)=It;
X(4,:,:)=Gt;
X(5,:,:)=HWt;
X(6,:,:)=IDP;
X(7,:,:)=RCt;
X(8,:,:)=ITAt;
X(9,:,:)=ICt;
X(10,:,:)=IAt;
X(11,:,:)=IRt;
X(12,:,:)=Rt;
X(13,:,:)=CRIt;
X(14,:,:)=TARIt;
X(15,:,:)=ARIt;
end

