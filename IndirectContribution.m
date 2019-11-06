function Ct = IndirectContribution(Yt,beta,Xt,indx,tau,EOCV,maxtau)
% Recursively computes the indirect contribution of the components sent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yt - the estimated attack rate
% beta - the coefficients for the regression model
% Xt - covariate
% indx - the index to compute the contribution
% tau - lags used in the regression model
% EOCV - the effectiveness of vaccination
% maxtau - the max lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ct - the temporal indirect effect of the compnenets in indx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ct=zeros(length(Yt(:,1)),length(Yt(1,:))+maxtau);

for ii=1:length(Yt(1,:))
    for jj=1:length(beta) % We do not include the first component which is just a constant
        if(~isempty(find(jj==indx)))
            Ct(:,ii+maxtau)=Ct(:,ii+maxtau)+(1-EOCV(:,ii)).*beta(jj).*(squeeze(Xt(jj,:,ii))'); 
        elseif(jj>1)
            Ct(:,ii+maxtau)=Ct(:,ii+maxtau)+(1-EOCV(:,ii)).*beta(jj).*(squeeze(Xt(jj,:,ii))').*Ct(:,ii+maxtau-tau(jj)); 
        end
    end
    Ct(:,ii+maxtau)=Ct(:,ii+maxtau)./Yt(:,ii);
end
Ct=Ct(:,(1+maxtau):end); % truncate to the points computed in the regression model

end

