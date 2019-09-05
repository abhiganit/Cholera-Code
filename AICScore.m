function AIC = AICScore(k,n,RSS)
% Returns the AIC score based on the least squares fit of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k - the number of paramters being estimated
% n - the number of data points used in the estimation
% RSS- Residual sum of squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIC - the AIC score

AIC=2*k+n*log(RSS);
end

