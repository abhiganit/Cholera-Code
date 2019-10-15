function Y = Et(WI,PS,IDPt)
% Computes the external incidence based on the wieghting matrix IDPt matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WI - Weekly incidence for the gov over time (Number Gov X NW)
% PS - Populatino size matrix for the outbreak (Number Gov X NW)
% IDPt - Matrix of the temporal change in IDP due ot conflict (size NW x
% Number Gov x Number Gov
    % IDPt(ii,jj,kk) - is the matrix for week ii, for gov jj where the
    % population went to from gov kk (i.e. people from kk went to jj on
    % week ii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Y - the external incidence (Number Gov X NW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
Y=zeros(size(WI)); % Allocate space for the output
NW=length(WI(1,:)); % number of weeks to compute
NG=length(WI(:,1)); % number of governorates


for ii=1:NW
    PIDP=squeeze(IDPt(ii,:,:));
    for jj=1:NG
        PIDP(:,jj)=PIDP(:,jj)./PS(jj,ii); % Dtermine the percentage of the population that displaced for the week (the column jj is the people leaving gov jj to gov ii)
    end
    Y(:,ii)=PIDP*WI(:,ii); % the external contribution
end


end

