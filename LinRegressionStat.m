function [SE,tStat,pValue] = LinRegressionStat(be,resid,DF,X,twotail)
% LinRegressionStat(be,resid,DF,X,twotail) produces the standard error
% t-staistic and pvalue for the beta coefficients we obtained from the
% fitting process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% be- the estimated beta value (size 1xM)
% resid- the vector of residual values
% DF - the number of parameters that were estimated during the fitting
% process
% X - The vector of X values for the beta value (NXM)
% twotail- Specify whether on wants a single tail p-value (two-tail=0) or the two
% tail (twotail=1) (vector size 1xM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SE - the Standard error
% tStat - the T- statistic
% pValue - The one sided p-value

%% Run computation
V=(sum(resid.^2)./(length(X(:,1))-DF)).*inv(X'*X); % The covariance matrix for the estimates
SE=sqrt(diag(V))'; % the standard error for the coefficients
tStat=be./SE; % The t-statistic for eahc coefficient
pValue=1-tcdf(tStat,length(X(:,1))-DF); % One-tail t-test value
f=find(twotail==1); % Find where you want to the two-tail test done
pValue(f)=pValue(f)+tcdf(-tStat(f),length(X(:,1))-DF); % Multiple by two since the distrituon is symettric

end

