function [F]= OFuncProGA(x,WI,tA,Ctv,Rtv,XU,maxtau,P,RC,H,WPIN,FPIN,Mt,Wheatt,Dieselt,V1,V2)
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
 % or one (1 x 18)
    % WPIN and FPIN - People in Need of WaSh and Food Security
    % PDG - population density for the govenrate
    % HFG- Health facilities in the govnorates
    % It - Incidence last week
    % Gt - Inicedence in the other govnorates
    % HWt - Health care and people in need
    % RCt - Rebel control
    % ITAt - Product of incidence (Indicator) and targeted attacks 
    % ICt - Product of incidence (Indicator) and conflict 
    % IAt - Product of incidence (Indicator) and air and shelling attacks 
    % IRt - Rainfall and wash (Indicator)
    % Rt- Rainfall incidence and WASH 
    % CRIt - Conflict, rianfall incidence, and WASH (Indicator)
    % TARIt -Targeted Attacks, rianfall, incidence and WASH (Indicator)
    % ARIt -Air/Shelling Attacks, rianfall, incidence and WASH (Indicator)
    % Wheatt - proce of wheat and food insecurity (Indicator)
    % Dieselt - proce of diesel and food insecurity (Indicator)
    % WaSHDieselt- price of diesel and WASH (Indicator)
% maxtau - maximum lag allowed for the model (used in the truncation of the
% data set
% P - population density for the governorates
% RC - rebel control
% H - the density of health facililities in the govnorates
% WPIN - wash people in need
% FPIN - food security people in need
% Mt - Sheeling and Ari attacks
% Wheatt - proce of wheat
% Dieselt- price of diesel
% V1 - proportion of the populatino to receive at least one dose
% V2 - proportion of the populatino to receive at least two doses
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
[xf] = ExpandPar(x,XU,0);
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF,mln,a,KV,dV]=RetParameterGA(xf,XU);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,~]= LogisticModel(beta,WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,FPIN,Mt,Wheatt,Dieselt,mln,a,V1,V2,KV,dV);

FF=(WI(:,(maxtau+1):end))-(Yt); % Compute the difference for the times and the locations that is tau weeks ahead
F=mean(sum(FF(:).^2)); % convert the matrix into a vector for the use of lsqnonlin
end

