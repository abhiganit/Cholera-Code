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
XU=[0 0 1 1 1 1]; 

% Specify the lag
% tau(1) - Past incidence
% tau(2) - Product of incidence and attacks
% tau(3) - Product of incidence and conflict
% tau(4) - Product of incidence and rainfall
% tau(5) - Perciptiation only

tau=[1 1 2 3 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The formation of the environmental function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify the attack function to be used
AF=1; % AF=0 attack only has effect before; AF=1 Attack has effect only after; AF=2; Attack has effect before and after

%Specify the conflict function to be used
CF=1; % CF=0 linear effect; CF=1 Hill function with n=1; CF=2; Full hill function; 

% Specify the rainfall function to be used
RF=2; % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Specify plots generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the govenorates want to examine (1-22)
%1- 'Abyan' 
%2- Aden' 
%3- Al Bayda' 
%4- Al Dhale''e' 
%5- Al Hudaydah' 
%6- Al Jawf' 
%7- Al Maharah' 
%8- Al Mahwit' 
%9- Amanat Al Asimah' 
%10- Amran' 
%11- Dhamar' 
%12- Hadramaut' 
%13- Hajjah' 
%14- Ibb' 
%15- Lahj' 
%16- Marib' 
%17- Raymah' 
%18- Sa''ada' 
%19- Sana''a' 
%20- Shabwah' 
%21- Socotra' 
%22- Taizz'
G=[]; 

% Specify if want plot of fit
PF=0; %PF=0 does not produce plot

% Specify if want plot of function
PE=1; %PE=0 does not produce plot

% Specify if want plot of projection and fit
PP=1; %PP=0 does not produce plot

% Specigy if table and fitting information displayed 
DT=1; % DT=0; does not display this information

% Run the fitting
RunFitting(XU,tau,AF,CF,RF,G,PF,PE,PP,DT);