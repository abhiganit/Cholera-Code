function [Yt,PDG,HFG,It,IAt,ICt,IRt,Rt,Gt,At,RCt]= LogisticModel(beta,WI,tA,DB,DA,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks
%===============================
% Input
%===============================
% t - the single time of interest
% It - The incidence for the different areas for the given week
% beta- the coefficients for the regression model
% tA- A  vector indicating the time of the attacks 
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
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
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form
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

%% Input for regression model

PDG=repmat(P,1,length(WI(1,(1+maxtau-tau(1)):(end-tau(1))))).*WI(:,(1+maxtau-tau(1)):(end-tau(1)));
HFG=repmat(H,1,length(WI(1,(1+maxtau-tau(2)):(end-tau(2))))).*WI(:,(1+maxtau-tau(2)):(end-tau(2)));
It=WI(:,(1+maxtau-tau(3)):(end-tau(3))); % Using the past incidence with a lag of tau weeks
IAt=WI(:,(1+maxtau-tau(4)):(end-tau(4))).*ImpactAttack(tA,DB,DA,tau(4),maxtau); % Product of incidence and attacks 
ICt=WI(:,(1+maxtau-tau(5)):(end-tau(5))).*ImpactConflict(Ctv(:,(1+maxtau-tau(5)):(end-tau(5))),K,n,CF); %Product of incidence and conflict
IRt=WI(:,(1+maxtau-tau(6)):(end-tau(6))).*ImpactRainfall(Rtv(:,(1+maxtau-tau(6)):(end-tau(6))),RIF,rl); %Product of incidence and rainfall
Rt=ImpactRainfall(Rtv(:,(1+maxtau-tau(7)):(end-tau(7))),RF,rh); %rainfall
Gt=repmat(sum(WI(:,(1+maxtau-tau(8)):(end-tau(8))),1),length(It(:,1)),1)-WI(:,(1+maxtau-tau(8)):(end-tau(8))); % Residual incidence (i.e. incdeicne in other govnorates)
At=ImpactAttack(tA,0,DAE,tau(9),maxtau);
RCt=repmat(RC,1,length(WI(1,(1+maxtau-tau(10)):(end-tau(10))))).*WI(:,(1+maxtau-tau(10)):(end-tau(10)));

%% Output of regression model: the predicted weekly incidence of the model
Yt=beta(1)+beta(2).*PDG+beta(3).*HFG+beta(4).*It+beta(5).*IAt+beta(6).*ICt+beta(7).*IRt+beta(8).*Rt+beta(9).*Gt+beta(10).*At+beta(11).*RCt; 

end

