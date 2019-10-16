function [k,beta,tau,DB,DA,DBE,DAE,K,n,rh,CF,RF]=RetParameterGA(x,XU)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k- the numebr of estimated paramters in the model
% beta- the coefficients for the regression model
% tau - the lag  (1X9)
% DB - the number of days before the attack that incidence is affected
% DA - the number of days after the attack that incidence is affected
% DBE - the number of days before the attack 
% DAE - the number of days after the attack that environment impacted
% K- the saturation function
% n - the hill coefficient
% rh - Threshold for rainfall for the covariate of rainfall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the coefficients fo rthe regression model

beta=[10.^x(1:length(XU))].*XU;
nob=length(XU);
tau=[1 x(nob+[1:8])];
CF=x(nob+[9:10]);
RF=x(nob+[11:14]);
lenbeta=length(XU)+14;
k=sum(XU)+ sum(XU([2:9])); % Count the number of coefficients being estimated the second sum is for estimating the lag of the different components


%% Attack asscoaited paramters
DA=zeros(2,1);
DB=zeros(1,1);

if(XU(3)==1)  
        DB(1)=10.^x(lenbeta+1); %looking befroe the attack
        k=k+1; % Add estimated paramter
%     end
        DA(1)=10.^x(lenbeta+2);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end

if(XU(8)==1)  % See if attacks being used at all

        DA(2)=10.^x(lenbeta+3);  %looking after the attack
        k=k+1; % Add estimated paramter

end

DAE=zeros(2,1);
DBE=zeros(1,1);

if(XU(5)==1)  % See if attacks being used at all

        DBE(1)=10.^x(lenbeta+4); %looking befroe the attack
        k=k+1; % Add estimated paramter

        DAE(1)=10.^x(lenbeta+5);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end

if(XU(9)==1)  % See if attacks being used at all

        DAE(2)=10.^x(lenbeta+6);  %looking after the attack
        k=k+1; % Add estimated paramter
%     end
end


%% Conflict associated paramters
K=zeros(2,1);
n=zeros(2,1);
if(XU(4)==1) % See if conflict is being used at alls
    if(CF(1)~=0)
        K(1)=10.^x(lenbeta+7); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(1)=0;
    end
    if(CF(1)==2) % If the full hill function is being used
        n(1)=10.^x(lenbeta+8); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(1)=1; % Hill coefficient ot estimated
    end    
end

if(XU(7)==1) % See if conflict is being used at alls
    if(CF(2)~=0)
        K(2)=10.^x(lenbeta+9); % Set rate of change for the paramter of the effects of conflict
        k=k+1; % add a paramter
    else
        K(2)=0;
    end
    if(CF(2)==2) % If the full hill function is being used
        n(2)=10.^x(lenbeta+10); % Hill coefficient estimate
        k=k+1; % add to estimated paramters
    else
        n(2)=1; % Hill coefficient ot estimated
    end    
end


 %% Rainfall assocaited paramters 
rh=zeros(4,1);
if(XU(6)==1) % See if rainfall is being used at all
    rh(1)=10.^x(lenbeta+11);
    k=k+1; % add paramrter    
end
if(XU(7)==1) % See if rainfall is being used at all
    rh(2)=10.^x(lenbeta+12);
    k=k+1; % add paramrter    
end
if(XU(8)==1) % See if rainfall is being used at all
    rh(3)=10.^x(lenbeta+13);
    k=k+1; % add paramrter    
end
if(XU(9)==1) % See if rainfall is being used at all
    rh(4)=10.^x(lenbeta+14);
    k=k+1; % add paramrter    
end

end

