function [k,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,DAR,w,sigma_w]=RetParameter_Table(x,XU,maxtau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - the log-10 parameters from the fitting process
% XU - Specify the model being fit    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k- the numebr of estimated paramters in the model
% beta- the coefficients for the regression model
% tau - the lag  
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
% DB - the number of days before the targeted attack that incidence is affected
% DA - the number of days after the targeted attack that incidence is affected
% DBE - the number of days before the air/shelling attack 
% DAE - the number of days after the air/shelling attack that environment impacted
% K- the saturation function
% n - the hill coefficient
% rl - Threshold for rainfall for the covariate of rainfall*indicator
% rh - Threshold for rainfall for the covariate of rainfall * incidence
% RIF - Rainfall function for rainfall* indicator
% mln - the threshold price for wheat, diesel
% a- the weight for the number of people in need WASH vs Food security
% KV - saturation value for the hill finction
% dV - deca yrate of the effect of vaccination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the coefficients fo rthe regression model

beta=[10.^x(1:length(XU))].*XU;
tau=repmat([1:maxtau],1,length(XU)/maxtau);
lenbeta=length(XU);
k=sum(XU); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
DA=NaN.*zeros(1,1);
DB=NaN.*zeros(1,1);
if(sum(XU(1:maxtau))>=1)  % See if attacks being used at all
    DA(1)=10.^x(lenbeta+1);  %looking after the attack
    k=k+1; % Add estimated paramter
end



%% Conflict associated paramters
K=NaN.*ones(2,1);
n=NaN.*ones(2,1);
if(sum(XU((maxtau+1):2*maxtau))>=1) % See if conflict is being used at alls
    
    K(1)=10.^x(lenbeta+2); % Set rate of change for the paramter of the effects of conflict
    k=k+1; % add a paramter

    n(1)=round(10.^x(lenbeta+3)); % Hill coefficient estimate
    k=k+1; % add to estimated paramters
end

% Shellings
if(sum(XU((2.*maxtau+1):3*maxtau))>=1)  % See if attacks being used at all
    
    K(2)=10.^x(lenbeta+4); % Set rate of change for the paramter of the effects of conflict
    k=k+1; % add a paramter

    n(2)=round(10.^x(lenbeta+5)); % Hill coefficient estimate
    k=k+1; % add to estimated paramters
end

KP=NaN.*zeros(2,1);
% Diesel price
if(sum(XU((3.*maxtau+1):4*maxtau))>=1)
    KP(1)=10.^x(lenbeta+6);
    k=k+1;
end

% Wheat price
if(sum(XU((4.*maxtau+1):5*maxtau))>=1)
    KP(2)=10.^x(lenbeta+7);
    k=k+1;
end

%Rainfall
r0=10.^x(lenbeta+8);
k=k+1;

%Rainfall
temp_0=10.^x(lenbeta+9);
k=k+1;

% Vaccination

KV=10.^x(lenbeta+[10]);
dV=[10.^x(lenbeta+11); exp(log(26/56)/(3*52)) ];
k=k+2;
w=10.^x(lenbeta+[12]);
k=k+1;
sigma_w=10.^x(lenbeta+[13]);
k=k+1;
DAR=10; % https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-018-6299-3 (50% asymp. 0.05-0.1 disease burden, then need to sautrate ontop by multiplying by 100 i.e. if all the population is infected the saturation function is less than 1%. Thus, if 50% cases are asymptomatic this implies 2 cases for every confirmed. Of the reported cases up to 0.05% are asscoaited wit hthe disease burded. thus, 0.1*100=10)
end