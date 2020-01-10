function [F]= MOFuncProGA(x,WI,tA,Ctv,Rtv,XU,maxtau,P,RC,H,WPIN,Mt)
% The difference of the predicted incidence and weekly incidence for all
% weeks and areas
%===============================
% Input
%===============================
% x - vector for the paramters being estimated in the model
% WI - Weekly incidence for the M areas for the N data points(MxN)
% tA- A vector indicating the time of the attacks
% Ctv - the vector for the conflict
% Rtv - the vector rainfall at time t
 % XU - the X input that will be used in the fitting process. Will be zeros
 % or one (1 x 11)
    % XU(1)- beta_0
    % XU(2) - population density
    % XU(3) - number of health facilities 
    % XU(4) - Past incidence
    % XU(5) - Product of incidence and attacks
    % XU(6) - Product of incidence and conflict
    % XU(7) - Product of cumulative attacks, incidence and rainfall
    % XU(8) - cumulative attacks and Rainfall        
    % XU(9) - Incidence in other govnorates
    % XU(10) - Attacks only
    % XU(11) - Rebel control
    % XU(12) - WASH
    % XU(13) - WASH and incidence
% tau - the lag
    % tau(1) - population density incidence
    % tau(2) - health zone incidence
    % tau(3) - Past incidence
    % tau(4) - Product of incidence and attacks
    % tau(5) - Product of incidence and conflict
    % tau(6) - Product of cumulative attacks, incidence and rainfall
    % tau(7) - cumulative attacks and Perciptiation
    % tau(8) - Incidence in other govneroates
    % tau(9)- Attack only
    % tau(10)- Rebel control
    % tau(11) - WASH and Incidence
% maxtau - maximum lag allowed for the model (used in the truncation of the
% data set
% AF - Specify the attack function to be used
    % AF=0 attack only has effect before; 
    % AF=1 Attack has effect only after; 
    % AF=2; Attack has effect before and after
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
% RIF - the indication of what function will be used for incidence*rainfall
    % RIF=0 the increase in incidence when rainfall is low
    % RIF=1 the increase in incidence when rainfall is high
    % RIF=2 the increase in incidence when rainfall is low and high
% RF - the indication of what function will be used for rainfall
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
% P - population density for the governorates
% RC - rebel control
% H - the density of health facililities in the govnorates
% WPIN - wash people in need
%=================================
% Output
%=================================
% F - the objective function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the paramter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find entries that are non-zero and non-one and set them to one

f=find(XU~=0); % find non-zero entries
g=find(XU(f)~=1); % among non-zero entries find the ones that are not one
XU(f(g))=1; % set non-zero and non-one to one


%Returns the paramters for the specified functions based on the
%transformation from the bounds
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(x,XU);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,~]= LogisticModel(beta,WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,Mt);

F=log10(sum(((WI(:,(maxtau+1):end))-(Yt)).^2,2)); % Compute the difference for the times and the locations that is tau weeks ahead
test=0;
end

