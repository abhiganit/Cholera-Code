function [k,beta,DB,DA,K,n,rl,rh]=RetParameter(x,XU,AF,CF,RF)
%Based on the input-x and functions used we return the proper paramters to
%evalaute the regression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x - the log-10 parameters from the fitting process
% XU - Specify the model being fit
    % XU(1)- beta_0
    % XU(2) - Past incidence
    % XU(3) - Product of incidence and attacks
    % XU(4) - Product of incidence and conflict
    % XU(5) - Product of incidence and rainfall
    % XU(6) - Rainfall only
    % XU(7) - Incidence in other govnorates
% AF - Specify the attack function to be used
    % AF=0 attack only has effect before; 
    % AF=1 Attack has effect only after; 
    % AF=2; Attack has effect before and after
% CF - what conflict function that is being used
        % CF=0 linear effect; 
        %CF=1 Hill function with n=1; 
        %CF=2; Full hill function
% RF - the indication of what function will be used
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k- the numebr of estimated paramters in the model
% beta- the coefficients for the regression model
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% K- the saturation function
% n - the hill coefficient
% rl - Rate the imapct in icidence decays from low rain fall
% rh - Rate the imapct in icidence decays from high rain fall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the coefficients fo rthe regression model
beta=[x(1) 10.^x(2:7)'].*XU;

k=sum(XU); % Count the number of coefficients being estimated

%% Attack asscoaited paramters
if(XU(3)>=1)  % See if attacks being used at all
    if(AF==1) % If only look after attack set DB=0
        DB=0; 
    else    
        DB=10.^x(8); %looking befroe the attack
        k=k+1; % Add estimated paramter
    end
    if(AF==0) % % If only look before attack set DA=0
        DA=0;
    else
        DA=10.^x(9);  %looking after the attack
        k=k+1; % Add estimated paramter
    end
else % Attack is not being used in the model
    DA=0;
    DB=0;
end


%% Conflict associated paramters
if(XU(4)>=1) % See if conflict is being used at alls
    K=10.^x(10); % Set rate of change for the paramter of the effects of conflict
    k=k+1; % add a paramter
    if(CF==2) % If the full hill function is being used
        n=10.^x(11); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n=1; % Hill coefficient ot estimated
    end    
else % Conflict is not being used in the model
    K=1;
    n=1;
end

 %% Rainfall assocaited paramters 
if(XU(5)>=1) % See if rainfall is being used at all
    rl=10.^x(12); % Esitmate for the low-rainfall heighted incidence
    if(RF~=1) % See if low-rainfall heighted incidence is being used
        k=k+1; % add paramrter
    end
    rh=10.^x(13); % Esitmate for the high-rainfall heighted incidence
    if(RF~=0) % See if high-rainfall heighted incidence is being used
        k=k+1;% add paramrter
    end
    
else % Rainfall is not being used at all
    rl=0;
    rh=0;
end

end

