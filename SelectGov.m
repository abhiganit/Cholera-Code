function [GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,PDS)

NGS=floor(length(GNZI)*PDS);
GTF=zeros(NGS,1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle

     % Include all places where vaccination occurs
     f=find(GV(GNZI)==1);
     GTF(1:sum(GV))=f;
     % Add remaining gov
    NGInc=sum(GTF>0);
    NumGovR=NGS-NGInc;
    if(NumGovR==1) % there are two fewer non-rebel vaccinated areas
        Itemp=sum(WI(GNZI,:),2);
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination

        % Find the top max
        for ii=(NGInc+1):(NGInc+ceil(NumGovR/4))
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
    elseif(NumGovR==2)
        Itemp=sum(WI(GNZI,:),2);
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination

        % Find the top max
        for ii=(NGInc+1):(NGInc+ceil(NumGovR/2))
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
        Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
        % Find the minimum contributors
        for ii=(NGInc+ceil(NumGovR/2)+1):(NGInc+ceil(NumGovR))
           f=find(Itemp==min(Itemp)); % Select minimum
           Itemp(f)=max(Itemp); % Set to maximum for not selected again
           GTF(ii)=f; % Record index
        end
    elseif(NumGovR==3)
        Itemp=sum(WI(GNZI,:),2);
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination

        % Find the top max
        for ii=(NGInc+1):(NGInc+1)
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
        Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
        % Find the minimum contributors
        for ii=(NGInc+2):(NGInc+2)
           f=find(Itemp==min(Itemp)); % Select minimum
           Itemp(f)=max(Itemp); % Set to maximum for not selected again
           GTF(ii)=f; % Record index
        end
        
        % Add Rebel control gov
        Itemp=sum(WI(GNZI,:),2);
        Itemp((RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination

        % Find the top max
        for ii=(NGInc+3):(NGInc+3)
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
    else
        % No Rebel control and no vaccination

        Itemp=sum(WI(GNZI,:),2);
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination

        % Find the top max
        for ii=(NGInc+1):(NGInc+ceil(NumGovR/4))
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
        Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
        Itemp((1-RC(GNZI)).*(1-GV(GNZI))==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
        % Find the minimum contributors
        for ii=(NGInc+ceil(NumGovR/4)+1):(NGInc+ceil(NumGovR/2))
           f=find(Itemp==min(Itemp)); % Select minimum
           Itemp(f)=max(Itemp); % Set to maximum for not selected again
           GTF(ii)=f; % Record index
        end

        % Rebel control and no vaccination

        Itemp=sum(WI(GNZI,:),2);
        Itemp((RC(GNZI)).*(1-GV(GNZI))==0)=0; % find the index where there is no rebel control or no vaccination vaccination
        % Find the top max
        for ii=(NGInc+ceil(NumGovR/2)+1):(NGInc+ceil(3.*NumGovR/4))
           f=find(Itemp==max(Itemp)); % Find the maximum
           Itemp(f)=0; % set maximum to zero for it is no longer selected
           GTF(ii)=f; % Record index
        end
        Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
        Itemp((RC(GNZI)).*(1-GV(GNZI))==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
        % Find the minimum contributors
        for ii=(NGInc+ceil(3.*NumGovR/4)+1):(NGInc+ceil(NumGovR))
           f=find(Itemp==min(Itemp)); % Select minimum
           Itemp(f)=max(Itemp); % Set to maximum for not selected again
           GTF(ii)=f; % Record index
        end
    end
GTF=sort(unique(GTF))'; % Gov. to used in the fitting of the model. We sort to keep order consistent with GNZI
GTCV=zeros(length(GNZI)-NGS,1); 
cc=1;
for ii=1:length(GNZI)
    if(isempty(find(GTF==ii)))
        GTCV(cc)=ii;
        cc=cc+1;
    end
end

end

