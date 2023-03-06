function [F]= OFuncProGAGLevel(x,WI,tA,Ctv,XU,maxtau,WPIN,FPIN,Mt,Wheatt,Dieselt,V1,V2,Rtv,Temptv,Pop,CI)
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
 % or one
% maxtau - maximum lag allowed for the model (used in the truncation of the
% data set
% P - population density for the governorates
% RC - rebel control
% H - the density of health facililities in the govnorates
% WPIN - wash people in need
% FPIN - food security people in need
% Mt - Sheeling and Ari attacks
% Wheatt - proce of wheat
% Dieselt- price of diesel
% V1 - proportion of the populatino to receive at least one dose
% V2 - proportion of the populatino to receive at least two doses
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
[xf] = ExpandPar(x,XU,maxtau);
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,temp_0,DAR,w,sigma_w]=RetParameterGA(xf,XU,maxtau);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,~]= LogisticModel(beta,tA,DB,DA,Ctv,K,n,tau,maxtau,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,V1,V2,KV,dV,Rtv,Temptv,r0,temp_0,WI,Pop,CI,DAR,w);

FF=log(normpdf(WI(:,(maxtau+1):end),Yt,sigma_w));
% FF=(WI(:,(maxtau+1):end))-(Yt); % Compute the difference for the times and the locations that is tau weeks ahead
F=-sum(FF(:)); % convert the matrix into a vector for the use of lsqnonlin
end

