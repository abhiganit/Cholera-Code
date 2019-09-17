function Y = ImpactAttack(tA,DB,DA,tau,maxtau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y = ImpactAttack(t,tA,DB,W,DE,NR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the impact of each attack based on the time of the outbreak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% t - the time to evaluate the function
% tA- A vector indicating the time of the attacks
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% tau- the lag being used in the regression model for all variables (Attack
% is the in location 2)
% maxtau- the maximum lag allowed to be used in the model such that all
% models use the same amount of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Y- the compunded effected of the attacks at time t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=tA; % day of the attack has full effect
MM=max(tA(:));
% See if we are examing days before or not
if(DB>0)
    WB=ceil((log(10^(-3))-log(MM))/log(DB)); % Look to see how many days vefore we need to go
    if(WB+1>length(tA(1,:))) % Ensure we do not exceed the length of the vector
       WB=length(tA(1,:))-1; % Set length based on the maximum otherwise
    end
else
    WB=0; % Not looking at the days before
end

% See if we are examing days after attack or not
if(DA>0)
    WA=ceil((log(10^(-3))-log(MM))/log(DA)); % Look to see how many days after we need to go
    if(WA+1>length(tA(1,:))) % Ensure we do not exceed the length of the vector
       WA=length(tA(1,:))-1; % Set length based on the maximum otherwise
    end
else
    WA=0; % Not looking at the days after
end

% Calcualte impact based on days before
for ii=1:WB % Number of days we need to go back
   temp=[tA(:,(ii+1):end) zeros(length(tA(:,1)),ii)]; % shift the matrix tA back to induce impact before the event occurrs
   Y=Y+(DB.^ii).*temp; % the level of impact decreases the further you go back as DB<1
end


% Calcualte impact based on days after the attack
for ii=1:WA
   temp=[zeros(length(tA(:,1)),ii) tA(:,1:(end-ii)) ]; % shift the matrix tA forward to induce impact after the event occurrs
   Y=Y+(DA.^ii).*temp; % the level of impact decreases the further ahead you go as DA<1
end


Y=Y(:,(1+maxtau-tau):(end-tau)); % truncate to integrate the lag of incidence
end

