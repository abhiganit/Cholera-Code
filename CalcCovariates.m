function [X] = CalcCovariates(WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN)
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
% X- the value of the covariates
% Constnat - ones 
% PDG - population density for the govenrate
% HFG- Health facilities in the govnorates
% It - Incidence last week
% IAt - Product of incidence and attacks 
% ICt - Product of incidence and conflict 
% IRt - Product of incidence and rainfall 
% Rt- Rainfall
% Gt - Inicedence in the other govnorates
% At- Attack only
% RCt - Rebel control

PDG=repmat(P,1,length(WI(1,(1+maxtau-tau(1)):(end-tau(1))))).*WI(:,(1+maxtau-tau(1)):(end-tau(1)));
HFG=repmat(H,1,length(WI(1,(1+maxtau-tau(2)):(end-tau(2))))).*WI(:,(1+maxtau-tau(2)):(end-tau(2)));
It=WI(:,(1+maxtau-tau(3)):(end-tau(3))); % Using the past incidence with a lag of tau weeks
IAt=WI(:,(1+maxtau-tau(4)):(end-tau(4))).*ImpactAttack(tA,DB,DA,tau(4),maxtau); % Product of incidence and attacks 
ICt=WI(:,(1+maxtau-tau(5)):(end-tau(5))).*ImpactConflict(Ctv(:,(1+maxtau-tau(5)):(end-tau(5))),K,n,CF); %Product of incidence and conflict
IRt=repmat(WPIN,1,length(WI(1,(1+maxtau-tau(6)):(end-tau(6))))).*WI(:,(1+maxtau-tau(6)):(end-tau(6))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(6)):(end-tau(6))),RIF,rl); %Product of incidence and rainfall
Rt=repmat(WPIN,1,length(WI(1,(1+maxtau-tau(7)):(end-tau(7))))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(7)):(end-tau(7))),RF,rh); %rainfall
Gt=repmat(sum(WI(:,(1+maxtau-tau(8)):(end-tau(8))),1),length(It(:,1)),1)-WI(:,(1+maxtau-tau(8)):(end-tau(8))); % Residual incidence (i.e. incdeicne in other govnorates)
At=ImpactAttack(tA,DBE,DAE,tau(9),maxtau);
RCt=repmat(RC,1,length(WI(1,(1+maxtau-tau(10)):(end-tau(10))))).*WI(:,(1+maxtau-tau(10)):(end-tau(10)));

X=zeros(11,length(PDG(:,1)),length(PDG(1,:)));

% Constant
X(1,:,:)=ones(size(PDG));
X(2,:,:)=PDG;
X(3,:,:)=HFG;
X(4,:,:)=It;
X(5,:,:)=IAt;
X(6,:,:)=ICt;
X(7,:,:)=IRt;
X(8,:,:)=Rt;
X(9,:,:)=Gt;
X(10,:,:)=At;
X(11,:,:)=RCt;
end

