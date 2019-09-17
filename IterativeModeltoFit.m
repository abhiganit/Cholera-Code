clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General model specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify XU the model being fit
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates


% Specify the lag
% tau(1) - Past incidence
% tau(2) - Product of incidence and attacks
% tau(3) - Product of incidence and conflict
% tau(4) - Product of incidence and rainfall
% tau(5) - Perciptiation only
% tau(6) - Incidence in other govneroates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The formation of the environmental function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify the attack function to be used
AF=1; % AF=0 attack only has effect before; AF=1 Attack has effect only after; AF=2; Attack has effect before and after

%Specify the conflict function to be used
CF=2; % CF=0 linear effect; CF=1 Hill function with n=1; CF=2; Full hill function; 

% Specify the rainfall function to be used
RF=2; % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall



% Run the iterative fitting
dW=4;
tau=[1 1 1 1 1 1];
XU=[1 0 0 1 0 0 0 0 0];
FittingIterative(XU,tau,AF,CF,RF,dW);
tau=[1 1 1 1 1 1];
XU=[0 1 0 1 0 0 0 0 0];
FittingIterative(XU,tau,AF,CF,RF,dW);
tau=[1 1 1 1 2 1];
XU=[1 0 0 1 0 0 0 1 0];
FittingIterative(XU,tau,AF,CF,RF,dW);
