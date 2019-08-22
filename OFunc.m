function F= OFunc(x,WI,A,C,R,tA,Ctv,Rtv,XU,tau,AF,CF,RF)
% The difference of the predicted incidence and weekly incidence for all
% weeks and areas
%===============================
% Input
%===============================
% x - vector for the paramters being estimated in the model
% WI - Weekly incidence for the M areas for the N data points(MxN)
% A- A=1 effect of attacks included A=0 No effects of attack
% C- C=1 effect of conflict included C=0 No effects of conflict
% R- R=1 effect of rainfall included R=0 No effects of rainfall
% tA- A vector indicating the time of the attacks
% Ctv - the vector for the conflict
% Rtv - the vector rainfall at time t
 % XU - the X input that will be used in the fitting process. Will be zeros
 % or one (1 x 6)
 % tau - the lag
% AF- the indication of what function is used for attacks
% CF - the indication of what function is used for conflict
% RF - the indication of what function will be used for rainfall
    % RF=0 the increase in incidence when rainfall is low
    % RF=1 the increase in incidence when rainfall is high
    % RF=2 the increase in incidence when rainfall is low and high
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
if(XU(5)==1)
    if(A~=0)
        A=1; % ensures that use of attacks either zero or one
    end
    if(C~=0)
        C=1; % ensures that use of conflict either zero or one
    end
    if(R~=0)
        R=1; % ensures that use of rainfall is either zero or one
    end
else
    A=0; % If no environmental factor then does not need to be included
    C=0; % If no environmental factor then does not need to be included
    R=0; % If no environmental factor then does not need to be included
end

%Returns the paramters for the specified functions based on the
%transformation from the bounds
[~,beta,DB,DA,K,n,rl,rh,ma,mc,mr]=RetParameter(x,XU,A,C,R,AF,CF,RF);

%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%
% Determine model predicted incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[Yt,~,~,~,~,~]=LogisticModel(beta,WI,A,tA,DB,DA,C,Ctv,K,n,R,Rtv,RF,rl,rh,tau,CF,ma,mc,mr); % Returns the incidence in matrix form of size Ng X (NW-tau)

F=WI(:,(max(tau)+1):end)-Yt; % Compute the difference for the times and the locations that is tau weeks ahead
F=F(:); % convert the matrix into a vector for the use of lsqnonlin
end

