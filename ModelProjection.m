function [Yt,Pt]= ModelProjection(beta,WI,tA,DB,DA,Ctv,K,n,Rtv,RF,rl,rh,tau,maxtau,CF,NWP)
% Produces the predicted incicence in matrix form for the diffrent areas
% and weeks and projects forward
%===============================
% Input
%===============================
% t - the single time of interest
% It - The incidence for the different areas for the given week
% beta- the coefficients for the regression model
% tA- A  vector indicating the time of the attacks 
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% Ct - the conflict at time t
% K- the saturation function
% n - the hill coefficient
% Rt - the rainfall at time t
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
% rl - Rate the imapct in icidence decays from low rain fall
% rh - Rate the imapct in icidence decays from high rain fall
% tau - the lag to use for the incidence and enviromental factors
    % tau(1) - Past incidence
    % tau(2) - Product of incidence and attacks
    % tau(3) - Product of incidence and conflict
    % tau(4) - Product of incidence and rainfall
    % tau(5) - Perciptiation only
    % tau(6) - Residual incidence
% maxtau- the maximum lag allowed to be used in the model such that all
% models use same information
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % NWP- the number of weeks to project forward
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form from the
% incidence data
% Pt- the projected incidence of the model in matrix form using model
% predicted incidence

%% Run the logistic model with the data
[Yt,~,~,~,~,~,~]= LogisticModel(beta,WI,tA(:,1:length(WI(1,:))),DB,DA,Ctv(:,1:length(WI(1,:))),K,n,Rtv(:,1:length(WI(1,:))),RF,rl,rh,tau,maxtau,CF);

%% Run the projection
temp=zeros(length(Yt(:,1)),1); % initialize the matrix for the projection
Pt=zeros(length(Yt(:,1)),NWP); % used for the 
for ii=1:NWP % Loop through the number of weeks that are to be projected
    WT=[WI temp]; % Need to append data to the end for the projection of incidence
    [temp2,~,~,~,~,~,~]= LogisticModel(beta,WT,tA(:,1:length(WT(1,:))),DB,DA,Ctv(:,1:length(WT(1,:))),K,n,Rtv(:,1:length(WT(1,:))),RF,rl,rh,tau,maxtau,CF); % Run model with appendend data
    Pt(:,ii)=temp2(:,end); % Record the projection of incidence
    temp=[Pt(:,1:ii) zeros(length(Yt(:,1)),1)];  % temporary variable to be appended
end

end

