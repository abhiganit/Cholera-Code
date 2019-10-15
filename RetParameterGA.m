function [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterGA(x,XU)
%Based on the input-x and functions used we return the proper paramters to
%evalaute the regression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - the log-10 parameters from the fitting process
% XU - Specify the model being fit    
        % XU(1)- beta_0
        % XU(2) - population density for the govenrate
        % XU(3)- Health facilities in the govnorates
        % XU(4) - Incidence last week
        % XU(5) - Inicedence in the other govnorates
        % XU(6) - Internally displaced people
        % XU(7) - Rebel control
        % XU(8) - Product of incidence and attacks 
        % XU(9) - Product of incidence and conflict 
        % XU(10) - Product of cumulative attacks incidence and rainfall 
        % XU(11)- cumulative attacks and Rainfall
        % XU(12) - Conflict, rianfall incidence
        % XU(13) - Attack, rianfall, incidence
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
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k- the numebr of estimated paramters in the model
% beta- the coefficients for the regression model
% tau - the lag  
    % tau(1) - Bias (not applicable)
    % tau(2) - population density incidence
    % tau(3) - health zone incidence
    % tau(4) - Past incidence
    % tau(5) - Incidence other governorates
    % tau(6) - Internally displaced
    % tau(7) - Rebel control
    % tau(8) - Product of incidence and attacks 
    % tau(9) - Product of incidence and conflict 
    % tau(10) - Product of cumulative attacks incidence and rainfall 
    % tau(11)- cumulative attacks and Rainfall
    % tau(12) - Conflict, rianfall incidence
    % tau(13) - Attack, rianfall, incidence
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% DBE - the number of days before the attack 
% DAE - the number of days after the attack that environment impacted
% K- the saturation function
% n - the hill coefficient
% rl - Threshold for rainfall for the covariate of rainfall*incidence
% rh - Threshold for rainfall for the covariate of rainfall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the coefficients fo rthe regression model

beta=[10.^x(1:length(XU))].*XU;
nob=length(XU);
tau=[1 1 1 1 1 1 1 x(nob+[1:8])];
CF=x(nob+[9:10]);
RF=x(nob+11);
RIF=x(nob+[12 13 14 15]);
lenbeta=length(XU)+15;
k=sum(XU)+ sum(XU([8:15])); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
DA=zeros(2,1);
DB=zeros(2,1);

if(XU(8)==1)  % See if attacks being used at all
%     if(AF==1) % If only look after attack set DB=0
%         DB=0; 
%     else    
        DB(1)=10.^x(lenbeta+1); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
%     if(AF==0) % % If only look before attack set DA=0
%         DA=0;
%     else
        DA(1)=10.^x(lenbeta+2);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end

if(XU(14)==1)  % See if attacks being used at all
%     if(AF==1) % If only look after attack set DB=0
%         DB=0; 
%     else    
        DB(2)=10.^x(lenbeta+3); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
%     if(AF==0) % % If only look before attack set DA=0
%         DA=0;
%     else
        DA(2)=10.^x(lenbeta+4);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end

DAE=zeros(2,1);
DBE=zeros(2,1);

if(XU(10)==1)  % See if attacks being used at all
%     if(AF==1) % If only look after attack set DB=0
%         DB=0; 
%     else    
        DBE(1)=10.^x(lenbeta+5); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
%     if(AF==0) % % If only look before attack set DA=0
%         DA=0;
%     else
        DAE(1)=10.^x(lenbeta+6);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end

if(XU(15)==1)  % See if attacks being used at all
%     if(AF==1) % If only look after attack set DB=0
%         DB=0; 
%     else    
        DBE(2)=10.^x(lenbeta+7); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
%     if(AF==0) % % If only look before attack set DA=0
%         DA=0;
%     else
        DAE(2)=10.^x(lenbeta+8);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
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
        K(2)=0;
    end
    if(CF(2)==2) % If the full hill function is being used
        n(2)=10.^x(lenbeta+12); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(2)=1; % Hill coefficient ot estimated
    end    
end


 %% Rainfall assocaited paramters 
rl=zeros(4,1);
if(XU(11)==1) % See if rainfall is being used at all
    rl(1)=10.^x(lenbeta+13);
    k=k+1; % add paramrter    
end
if(XU(13)==1) % See if rainfall is being used at all
    rl(2)=10.^x(lenbeta+14);
    k=k+1; % add paramrter    
end
if(XU(14)==1) % See if rainfall is being used at all
    rl(3)=10.^x(lenbeta+15);
    k=k+1; % add paramrter    
end
if(XU(15)==1) % See if rainfall is being used at all
    rl(4)=10.^x(lenbeta+16);
    k=k+1; % add paramrter    
end

 %% Rainfall only
if(XU(12)==1) % See if rainfall is being used at all
    rh=10.^x(lenbeta+17);
    k=k+1; % add paramrter
    
else % Rainfall is not being used at all
    rh=0;
end


end

