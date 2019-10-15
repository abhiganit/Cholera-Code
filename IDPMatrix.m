function IDPt = IDPMatrix(IDP,SNames,NW)
%IDPMATRIX constructs the temporal change in IDP matrix for the different
%areas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDP - Location dispplaced, Where came from, Year of discplacment, Month
% of displacement, number discplaced due to conflict
% SNames - the names of the govnerorates
% NW - the number of weeks that we want to go out to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDPt - Matrix of the temporal change in IDP due ot conflict (size NW x
% length(SNames) x length(SNames)
    % IDPt(ii,jj,kk) - is the matrix for week ii, for gov jj where the
    % population went to from gov kk (i.e. people from kk went to jj on
    % week ii)
StartIndx= datenum('10-03-2016');% The week of our first data point (October 3, 2016)
IDPt=zeros(NW,length(SNames),length(SNames)); % allocate space of the matirx
NN=length(SNames); % the number of areas examingn
for ii=1:NW % loop through the weeks    
   M=str2num(datestr(StartIndx+7*(ii-1),'mm')); % Month looking for
   Y=str2num(datestr(StartIndx+7*(ii-1),'yyyy')); % Year looking for
   if((Y==2019)) % We do not have data for after November 2018 so we will compute the monthly average 2019 year based on 2016-2018 (3 years)  2019 data goes up to Sept. So we divide by 3    
       fm=find(cell2mat(IDP(:,4))==M); % find the month
       fmfy=fm; % the rows containg the month for 2016-2018
       for wt=1:NN % Where the IDP went to
          fwt=find(strcmp(IDP(fmfy,1),SNames(wt))); % Search where they went to for gov wt
          fwtfmfy=fmfy(fwt); % update the index
          if(~isempty(fwtfmfy)) % only search if non-empty
              for cf=1:NN % Where the IDP came from
                fcf=find(strcmp(IDP(fwtfmfy,2),SNames(cf))); % search where the came from for gov cf
                fcffwtfmfy=fwtfmfy(fcf); % update index
                IDPt(ii,wt,cf)=sum(cell2mat(IDP(fcffwtfmfy,5)))./3; % number of people displaced due to conflict average for the three years
              end
          end
          IDPt(ii,wt,wt)=0; % Set the IDP to itself to zero as we are interested in using this to capture the spread of incidence among external areas
       end
   elseif((Y==2018)&&(M==12))  % We do not have data for after November 2018 so we will compute the monthly average 2019 year based on 2016-2017 (2 years)  
       fm=find(cell2mat(IDP(:,4))==M); % find the month
       fmfy=fm; % the rows containg the month for 2016-2018
       for wt=1:NN % Where the IDP went to
          fwt=find(strcmp(IDP(fmfy,1),SNames(wt))); % Search where they went to for gov wt
          fwtfmfy=fmfy(fwt); % update the index
          if(~isempty(fwtfmfy)) % only search if non-empty
              for cf=1:NN % Where the IDP came from
                fcf=find(strcmp(IDP(fwtfmfy,2),SNames(cf))); % search where the came from for gov cf
                fcffwtfmfy=fwtfmfy(fcf); % update index
                IDPt(ii,wt,cf)=sum(cell2mat(IDP(fcffwtfmfy,5)))./2; % number of people displaced due to conflict average for the three years
              end
          end
          IDPt(ii,wt,wt)=0; % Set the IDP to itself to zero as we are interested in using this to capture the spread of incidence among external areas
       end
   else % Anything before December 2018              
       fm=find(cell2mat(IDP(:,4))==M); % find the month
       fy=find(cell2mat(IDP(fm,3))==Y); % find the year
       fmfy=fm(fy); % the rows containg the month and year of interest
       for wt=1:NN % Where the IDP went to
          fwt=find(strcmp(IDP(fmfy,1),SNames(wt))); % Search where they went to for gov wt
          fwtfmfy=fmfy(fwt); % update the index
          if(~isempty(fwtfmfy)) % only search if non-empty
              for cf=1:NN % Where the IDP came from
                fcf=find(strcmp(IDP(fwtfmfy,2),SNames(cf))); % search where the came from for gov cf
                fcffwtfmfy=fwtfmfy(fcf); % update index
                IDPt(ii,wt,cf)=sum(cell2mat(IDP(fcffwtfmfy,5))); % number of people displaced due to conflict
              end
          end
          IDPt(ii,wt,wt)=0; % Set the IDP to itself to zero as we are interested in using this to capture the spread of incidence among external areas
       end

   end
end

end

