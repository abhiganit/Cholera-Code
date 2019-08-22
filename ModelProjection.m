function [Yt,Pt]= ModelProjection(beta,WI,A,tA,DB,DA,C,Ctv,K,n,R,Rtv,RF,rl,rh,tau,CF,ma,mc,mr,NWP)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks and projects forward
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
 % NWP- the number of weeks to project forward
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form from the
% incidence data
% Pt- the projected incidence of the model in matrix form using model
% predicted incidence

%% Run the logistic model with the data
[Yt,~,~,~,~,~]= LogisticModel(beta,WI,A,tA(:,1:length(WI(1,:))),DB,DA,C,Ctv(:,1:length(WI(1,:))),K,n,R,Rtv(:,1:length(WI(1,:))),RF,rl,rh,tau,CF,ma,mc,mr);

%% Run the projection
temp=zeros(length(Yt(:,1)),1); % initialize the matrix for the projection
Pt=zeros(length(Yt(:,1)),NWP); % used for the 
for ii=1:NWP % Loop through the number of weeks that are to be projected
    WT=[WI temp]; % Need to append data to the end for the projection of incidence
    [temp2,~,~,~,~,~]= LogisticModel(beta,WT,A,tA(:,1:length(WT(1,:))),DB,DA,C,Ctv(:,1:length(WT(1,:))),K,n,R,Rtv(:,1:length(WT(1,:))),RF,rl,rh,tau,CF,ma,mc,mr); % Run model with appendend data
    Pt(:,ii)=temp2(:,end); % Record the projection of incidence
    temp=[Pt(:,1:ii) zeros(length(Yt(:,1)),1)];  % temporary variable to be appended
end

end

