function F= OFunc(x,WI,tA,Ctv,Rtv,XU,tau,maxtau,AF,CF,RF,P,H)
% The difference of the predicted incidence and weekly incidence for all
% weeks and areas
%===============================
% Input
%===============================
% x - vector for the paramters being estimated in the model
% WI - Weekly incidence for the M areas for the N data points(MxN)
% tA- A vector indicating the time of the attacks
% Ctv - the vector for the conflict
% Rtv - the vector rainfall at time t
 % XU - the X input that will be used in the fitting process. Will be zeros
 % or one (1 x 6)
    % XU(1)- beta_0
    % XU(2) - Past incidence
    % XU(3) - Product of incidence and attacks
    % XU(4) - Product of incidence and conflict
    % XU(5) - Product of incidence and rainfall
    % XU(6) - Rainfall only
% tau - the lag
    % tau(1) - Past incidence
    % tau(2) - Product of incidence and attacks
    % tau(3) - Product of incidence and conflict
    % tau(4) - Product of incidence and rainfall
    % tau(5) - Perciptiation only
% maxtau - maximum lag allowed for the model (used in the truncation of the
% data set
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
% P - population density for the governorates
% H - the density of health facililities in the govnorates
%=================================
% Output
%=================================
% F - the objective function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the paramter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find entries that are non-zero and non-one and set them to one

f=find(XU~=0); % find non-zero entries
g=find(XU(f)~=1); % among non-zero entries find the ones that are not one
XU(f(g))=1; % set non-zero and non-one to one


%Returns the paramters for the specified functions based on the
%transformation from the bounds
[~,beta,DB,DA,K,n,rl,rh]=RetParameter(x,XU,AF,CF,RF);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,~,~,~,~,~,~,~]=LogisticModel(beta,WI,tA,DB,DA,Ctv,K,n,Rtv,RF,rl,rh,tau,maxtau,CF,P,H); % Returns the incidence in matrix form of size Ng X (NW-tau)

F=WI(:,(maxtau+1):end)-Yt; % Compute the difference for the times and the locations that is tau weeks ahead
F=F(:); % convert the matrix into a vector for the use of lsqnonlin
end

