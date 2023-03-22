function LL=log_likelihood_function(maxtau,WI,Yt,sigma_w)

% re-scale to purely a per-capita bases and not per 10,000
Data_Y=WI(:,(maxtau+1):end)./10000;
Model_Y=Yt./10000;

LL=log(normpdf(Data_Y,Model_Y,sigma_w));

% Data_Y=WI(:,(maxtau+1):end);
% Model_Y=Yt;
% LL=(Data_Y(:)-Model_Y(:)).^2;
% LL=-LL; % As objective function is set of for maximiziing log-likelihood
end