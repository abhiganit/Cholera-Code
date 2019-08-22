function Y = ImpactRainfall(Rt,RF,rl,rh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y = ImpactRainfall(Rt,H,rl,rh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the impact of rainfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Rt - the rainfall at time t
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
% rl - Rate the imapct in icidence decays from low rain fall
% rh - Rate the imapct in icidence decays from high rain fall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Y- Reutnrs the impact for the amount of rainfall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch RF
    case 0 % the increase in incidence when rainfall is low
        Y=exp(-rl.*Rt);
    case 1 % the increase in incidence when rainfall is high
        Y=exp(rh.*Rt)-1;
    case 2 % the increase in incidence when rainfall is high
        Y=exp(-rl.*Rt)+exp(rh.*Rt)-1;
    otherwise
        warning('Error: Appropriate function not selected');
        return;
end
Y(Y>1)=1; % Make sure that the measure does not exceed 100%
end
