function [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF,mln,a,KV,dV]=RetParameterPS(x,XU)
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
nob=length(XU);
tau=[ones(1,7) ceil(4.*x(nob+[1:(length(XU)-7)]))];
NTE=(length(XU)-7);
CF=ceil(3.*x(nob+NTE+[1 2]))-1;
RF=ceil(3.*x(nob+NTE+3))-1;
RIF=ceil(3.*x(nob+NTE+[4 5 6 7]))-1;
lenbeta=length(XU)+NTE+7;
k=sum(XU)+ sum(XU([8:length(XU)])); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
DA=zeros(2,1);
DB=zeros(2,1);
if(XU(8)==1)  % See if attacks being used at all
    DB(1)=10.^x(lenbeta+1); %looking befroe the attack
    k=k+1; % Add estimated paramter
    DA(1)=10.^x(lenbeta+2);  %looking after the attack
    k=k+1; % Add estimated paramter
end

if(XU(14)==1)  % See if attacks being used at all
        DB(2)=10.^x(lenbeta+3); %looking befroe the attack
        k=k+1; % Add estimated paramter

        DA(2)=10.^x(lenbeta+4);  %looking after the attack
        k=k+1; % Add estimated paramter
end

DAE=zeros(2,1);
DBE=zeros(2,1);

if(XU(10)==1)  % See if attacks being used at all
 
        DBE(1)=10.^x(lenbeta+5); %looking befroe the attack
        k=k+1; % Add estimated paramter
        DAE(1)=10.^x(lenbeta+6);  %looking after the attack
        k=k+1; % Add estimated paramter

end

if(XU(15)==1)  % See if attacks being used at all

        DBE(2)=10.^x(lenbeta+7); %looking befroe the attack
        k=k+1; % Add estimated paramter
        DAE(2)=10.^x(lenbeta+8);  %looking after the attack
        k=k+1; % Add estimated paramter
end



%% Conflict associated paramters
K=zeros(2,1);
n=zeros(2,1);
if(XU(9)==1) % See if conflict is being used at alls
    if(CF(1)~=0)
        K(1)=10.^x(lenbeta+9); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(1)=0;
    end
    if(CF(1)==2) % If the full hill function is being used
        n(1)=10.^x(lenbeta+10); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(1)=1; % Hill coefficient ot estimated
    end    
end

if(XU(13)==1) % See if conflict is being used at alls
    if(CF(2)~=0)
        K(2)=10.^x(lenbeta+11); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(1)=0;
    end
    if(CF(2)==2) % If the full hill function is being used
        n(1)=10.^x(lenbeta+12); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(2)=1; % Hill coefficient ot estimated
    end    
end


 %% Rainfall assocaited paramters 
rl=zeros(4,1);
if(XU(11)>=1) % See if rainfall is being used at all
    rl(1)=10.^x(lenbeta+13);
    k=k+1; % add paramrter    
end
if(XU(13)>=1) % See if rainfall is being used at all
    rl(2)=10.^x(lenbeta+14);
    k=k+1; % add paramrter    
end
if(XU(14)>=1) % See if rainfall is being used at all
    rl(3)=10.^x(lenbeta+15);
    k=k+1; % add paramrter    
end
if(XU(15)>=1) % See if rainfall is being used at all
    rl(4)=10.^x(lenbeta+16);
    k=k+1; % add paramrter    
end

 %% Rainfall only
if(XU(12)>=1) % See if rainfall is being used at all
    rh=10.^x(lenbeta+17);
    k=k+1; % add paramrter
    
else % Rainfall is not being used at all
    rh=0;
end

% Threshold for Wheat
mln=zeros(3,1);
if(XU(16)>=1)
    mln(1)=10.^x(lenbeta+18);
    k=k+1;
end
% Price threshold diesel
if(XU(17)>=1)
    mln(2)=10.^x(lenbeta+19);
    k=k+1;
end
% Price threshold diesel
if(XU(18)>=1)
    mln(3)=10.^x(lenbeta+20);
    k=k+1;
end

% Weight of WASH vs Food security
a=10.^x(lenbeta+21);
k=k+1;

% Vaccination
KV=10.^x(lenbeta+[22:23]);
dV=10.^x(lenbeta+24);
k=k+3;
end
