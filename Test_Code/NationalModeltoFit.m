clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General model specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify XU the model being fit
% XU(1)- beta_0
% XU(2) - Past incidence
% XU(3) - Product of incidence and attacks
% XU(4) - Product of incidence and conflict
% XU(5) - Product of incidence and rainfall
% XU(6) - Rainfall only
XU=[1 0 0 1 0 0]; 

% Specify the lag
% tau(1) - Past incidence
% tau(2) - Product of incidence and attacks
% tau(3) - Product of incidence and conflict
% tau(4) - Product of incidence and rainfall
% tau(5) - Perciptiation only

tau=[1 1 1 1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The formation of the environmental function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify the attack function to be used
AF=1; % AF=0 attack only has effect before; AF=1 Attack has effect only after; AF=2; Attack has effect before and after

%Specify the conflict function to be used
CF=0; % CF=0 linear effect; CF=1 Hill function with n=1; CF=2; Full hill function; 

% Specify the rainfall function to be used
RF=0; % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Specify plots generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specify if want plot of fit
PF=1; %PF=0 does not produce plot

% Specify if want plot of function
PE=1; %PE=0 does not produce plot


% Specigy if table and fitting information displayed 
DT=1; % DT=0; does not display this information

% Run the fitting

    for ii=0:1
        for ind=4:5
            for tt=1:4
                tau=[1 1 1 1 1];
                XU=zeros(1,6);
                XU(1)=ii;
                XU(3)=1;
                XU(6)=1;
                XU(ind)=1;
                tau(2)=1;
                tau(5)=2;
                tau(ind-1)=tt;
                RunFittingNational(XU,tau,AF,CF,RF,PF,PE,DT);
            end
        end
    end