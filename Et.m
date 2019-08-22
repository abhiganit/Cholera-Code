function Y = Et(A,tA,DB,DA,C,CF,Ct,K,n,R,Rt,H,rl,rh,ma,mc,mr,tau)
% Environment function that will be used in the regression modell
%===============================
% Input
%===============================
% t - the single time of interest
% A- A=1 effect of attacks included A=0 No effects of attack
% tA- A vector indicating the time of the attacks
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% C- C=1 effect of conflict included C=0 No effects of conflict
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
% Ct - the conflict at time t
% K- the saturation function
% n - the hill coefficient
% R- R=1 effect of rainfall included R=0 No effects of rainfall
% Rt - the rainfall at time t
% H - the indication of what function will be used
    % H=0 the increase in incidence when rainfall is low
    % H=1 the increase in incidence when rainfall is high
    % H=2 the increase in incidence when rainfall is low and high
% rl - Rate the imapct in icidence decays from low rain fall
% rh - Rate the imapct in icidence decays from high rain fall
 % ma - the magnitdue of attacks in the Envirmental factor compounded
 % effect
 % mc - the magnitdue of conflict in the Envirmental factor compounded
 % effect
 % mr - the magnitdue of rainfall in the Envirmental factor compounded
 % effect
 % tau - the lag being used in the model for all variables
%=================================
% Output
%=================================
% Y - The value of the environment based on the measures specified

% Initialize the impact of these different effects
Y=1;

if(A==1) % if attacks have an effect
   Y=Y.*(1-ma.*ImpactAttack(tA,DB,DA,tau)); % The compunded effects of the attacks specified
end
if(C==1) % if the daily level of conflict has an effect
    Y = Y.*(1-mc.*ImpactConflict(Ct,K,n,CF)); % The effects of the daily conflict levels specified
end
if(R==1)% if rainfall has an effect
    Y=Y.*(1-mr.*ImpactRainfall(Rt,H,rl,rh));
end  
Y=1-Y; % 1-product of the compunded effect 
end

