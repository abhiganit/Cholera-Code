function [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterGA(x,XU)
%Based on the input-x and functions used we return the proper paramters to
%evalaute the regression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - the log-10 parameters from the fitting process
% XU - Specify the model being fit    
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates
        % XU(10) - Attacks only
        % XU(11) - Rebel control
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
    % tau(11)- WASH status and incidence
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
tau=[1 1 1 x(nob+1) x(nob+2) x(nob+3) x(nob+4) 1 x(nob+5) 1 x(nob+6) x(nob+7)];
CF=x(nob+8);
RF=x(nob+9);
RIF=x(nob+10);
lenbeta=length(XU)+10;
k=sum(XU)+ sum(XU([5 6 7 8 10 12 13])); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
if(sum(XU([5 13]))>=1)  % See if attacks being used at all
%     if(AF==1) % If only look after attack set DB=0
%         DB=0; 
%     else    
        DB=10.^x(lenbeta+1); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
%     if(AF==0) % % If only look before attack set DA=0
%         DA=0;
%     else
        DA=10.^x(lenbeta+2);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
else % Attack is not being used in the model
    DA=0;
    DB=0;
end


%% Conflict associated paramters
if(sum(XU([6 12]))>=1) % See if conflict is being used at alls
    if(CF~=0)
        K=10.^x(lenbeta+3); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K=0;
    end
    if(CF==2) % If the full hill function is being used
        n=10.^x(lenbeta+4); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n=1; % Hill coefficient ot estimated
    end    
else % Conflict is not being used in the model
    K=1;
    n=1;
end

 %% Rainfall assocaited paramters 
if(XU(7)>=1) % See if rainfall is being used at all
    rl=10.^x(lenbeta+5);
    k=k+1; % add paramrter    
else % Rainfall is not being used at all
    rl=0;
end

 %% Rainfall only
if(XU(8)>=1) % See if rainfall is being used at all
    rh=10.^x(lenbeta+6);
    k=k+1; % add paramrter
    
else % Rainfall is not being used at all
    rh=0;
end

 %% Attack only
if(XU(10)>=1) % See if attack only being used
    DAE=10.^x(lenbeta+7);
    k=k+1; % add paramrter
    DBE=10.^x(lenbeta+8);
    k=k+1; % add paramrter
else % Rainfall is not being used at all
    DAE=0;
    DBE=0;
end

end

