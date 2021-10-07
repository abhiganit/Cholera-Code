function Ct = IndirectContribution(Yt,beta,Xt,indx,EOCV)
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

Ct=zeros(size(Yt));

for ii=1:length(Yt(1,:))
    for jj=1:length(indx) % We do not include the first component which is just a constant        
       if(ii>1)
        Ct(:,ii)=Ct(:,ii)+(sum((1-EOCV(:,1:ii)).*beta(indx(jj)).*(squeeze(Xt(indx(jj),:,1:ii))),2)./sum(Yt(:,1:ii),2));
       else
           Ct(:,ii)=Ct(:,ii)+(((1-EOCV(:,1:ii)).*beta(indx(jj)).*(squeeze(Xt(indx(jj),:,1:ii)))')./(Yt(:,1:ii)));
       end
    end
end

end

