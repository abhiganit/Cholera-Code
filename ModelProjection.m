function [Yt,Pt]= ModelProjection(beta,WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,NWP)
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
% DBE - the number of days before the attack that environemnt is affected
% DAE - the number of days after the attack that environemnt is affected
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
    % tau(1) - population density incidence
    % tau(2) - health zone incidence
    % tau(3) - Past incidence
    % tau(4) - Product of incidence and attacks
    % tau(5) - Product of incidence and conflict
    % tau(6) - Product of incidence and rainfall
    % tau(7) - Perciptiation only
    % tau(8) - Incidence in other govneroates
    % tau(9)- Attack only
    % tau(10)- Rebel control
% maxtau- the maximum lag allowed to be used in the model such that all
% models use same information
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
 % P - population density
 % RC - Reble control
 % H- health zones
 % NWP- the number of weeks to project forward
%=================================
% Output
%=================================
% Yt - the predicted incidence of the model in matrix form from the
% incidence data
% Pt- the projected incidence of the model in matrix form using model
% predicted incidence

%% Run the logistic model with the data
[Yt,~]= LogisticModel(beta,WI,tA(:,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(:,1:length(WI(1,:))),K,n,Rtv(:,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H);

%% Run the projection
temp=zeros(length(Yt(:,1)),1); % initialize the matrix for the projection
Pt=zeros(length(Yt(:,1)),NWP); % used for the 
for ii=1:NWP % Loop through the number of weeks that are to be projected
    WT=[WI temp]; % Need to append data to the end for the projection of incidence
    [temp2,~]= LogisticModel(beta,WT,tA(:,1:length(WT(1,:))),DB,DA,DBE,DAE,Ctv(:,1:length(WT(1,:))),K,n,Rtv(:,1:length(WT(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H); % Run model with appendend data
    Pt(:,ii)=temp2(:,end); % Record the projection of incidence
    temp=[Pt(:,1:ii) zeros(length(Yt(:,1)),1)];  % temporary variable to be appended
end

end

