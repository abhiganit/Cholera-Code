function [Yt,It,IAt,ICt,IRt,IEt]= LogisticModel(beta,WI,A,tA,DB,DA,C,Ctv,K,n,R,Rtv,RF,rl,rh,tau,CF,ma,mc,mr)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks
%===============================
% Input
%===============================
% t - the single time of interest
% It - The incidence for the different areas for the given week
% beta- the coefficients for the regression model
% A- A=1 effect of attacks included A=0 No effects of attack
% tA- A  vector indicating the time of the attacks 
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% C- C=1 effect of conflict included C=0 No effects of conflict
% Ct - the conflict at time t
% K- the saturation function
% n - the hill coefficient
% R- R=1 effect of rainfall included R=0 No effects of rainfall
% Rt - the rainfall at time t
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
% rl - Rate the imapct in icidence decays from low rain fall
% rh - Rate the imapct in icidence decays from high rain fall
% tau - the lag to use for the incidence and enviromental factors
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % ma - the magnitdue of attacks in the Envirmental factor compounded
 % effect
 % mc - the magnitdue of conflict in the Envirmental factor compounded
 % effect
 % mr - the magnitdue of rainfall in the Envirmental factor compounded
 % effect
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form
% It - Incidence last week
% IAt - Product of incidence and attacks 
% ICt - Product of incidence and conflict 
% IRt - Product of incidence and rainfall 
% IEt - Product of incidence and Invironmental factors

%% Input for regression model

It=WI(:,(1+max(tau)-tau(1)):(end-tau(1))); % Using the past incidence with a lag of tau weeks
IAt=WI(:,(1+max(tau)-tau(2)):(end-tau(2))).*ImpactAttack(tA,DB,DA,tau); % Product of incidence and attacks 
ICt=WI(:,(1+max(tau)-tau(3)):(end-tau(3))).*ImpactConflict(Ctv(:,(1+max(tau)-tau(3)):(end-tau(3))),K,n,CF); %Product of incidence and conflict
IRt=WI(:,(1+max(tau)-tau(4)):(end-tau(4))).*ImpactRainfall(Rtv(:,(1+max(tau)-tau(4)):(end-tau(4))),RF,rl,rh); %Product of incidence and rainfall
IEt=WI(:,(1+max(tau)-tau(5)):(end-tau(5))).*Et(A,tA,DB,DA,C,CF,Ctv(:,(1+max(tau)-tau(5)):(end-tau(5))),K,n,R,Rtv(:,(1+max(tau)-tau(5)):(end-tau(5))),RF,rl,rh,ma,mc,mr,[max(tau) tau(5)]); %Product of incidence and Invironmental factors

%% Output of regression model: the predicted weekly incidence of the model
Yt=beta(1)+beta(2).*It+beta(3).*IAt+beta(4).*ICt+beta(5).*IRt+beta(6).*IEt; 

end

