function [k,beta,tau,DB,DA,K,n,KP,KV,dV,r,r0,rm]=RetParameterGA(x,XU,CF,RF)
%Based on the input-x and functions used we return the proper paramters to
%evalaute the regression model

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
% CF - conflict function
% RIF - Rainfall function for rainfall* indicator
% RF - rainfall function for rainfall * incidence
% mln - the threshold price for wheat, diesel
% a- the weight for the number of people in need WASH vs Food security
% KV - saturation value for the hill finction
% dV - deca yrate of the effect of vaccination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the coefficients fo rthe regression model

beta=[10.^x(1:length(XU))].*XU;
tau=[ones(1,length(XU))];
lenbeta=length(XU);
k=sum(XU); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
DA=zeros(1,1);
DB=zeros(1,1);
if(XU(2)==1)  % See if attacks being used at all
    DA(1)=10.^x(lenbeta+1);  %looking after the attack
    k=k+1; % Add estimated paramter
end



%% Conflict associated paramters
K=zeros(4,1);
n=ones(4,1);
if(XU(3)==1) % See if conflict is being used at alls
    if(CF(1)==2)
        K(1)=10.^x(lenbeta+2); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(1)=0;
    end
    if(CF(1)~=0) % If the full hill function is being used
        n(1)=10.^x(lenbeta+3); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(1)=1; % Hill coefficient ot estimated
    end    
end

if(XU(7)==1) % See if conflict is being used at alls
    if(CF(2)==2)
        K(2)=10.^x(lenbeta+4); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(2)=0;
    end
    if(CF(2)~=0) % If the full hill function is being used
        n(2)=10.^x(lenbeta+5); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(2)=1; % Hill coefficient ot estimated
    end    
end
% Shellings
if(XU(4)==1)  % See if attacks being used at all
    if(CF(1)==2)
        K(3)=10.^x(lenbeta+6); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(3)=0;
    end
    if(CF(1)~=0) % If the full hill function is being used
        n(3)=10.^x(lenbeta+7); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(3)=1; % Hill coefficient ot estimated
    end
end

if(XU(8)==1)   % See if attacks being used at all
    if(CF(2)==2)
        K(4)=10.^x(lenbeta+8); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(4)=0;
    end
    if(CF(2)~=0) % If the full hill function is being used
        n(4)=10.^x(lenbeta+9); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(4)=1; % Hill coefficient ot estimated
    end
end

KP=zeros(2,1);
% Diesel price
if(sum(XU([5 9])>=1))
    KP(1)=10.^x(lenbeta+10);
    k=k+1;
end

% Wheit price
if(XU(10)>=1)
    KP(2)=10.^x(lenbeta+11);
    k=k+1;
end

% Vaccination

KV=10.^x(lenbeta+[12]);
dV=[10.^x(lenbeta+13); exp(log(26/56)/(4*52)) ];
k=k+2;

r=10.^x(lenbeta+[14]);
k=k+1;

if(RF(1)>=0)
    r0=10.^x(lenbeta+[15]);
    k=k+1;
else
    r0=0;
end
if(RF(2)>=0)
    rm=10.^x(lenbeta+[16]);
    k=k+1;
else
    rm=0;
end
end
