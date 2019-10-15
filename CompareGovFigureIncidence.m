function CompareGovFigureIncidence(Yt,Yt2,WID,GNZI,GTF)
% Plots the incidence data and model prediction for the area of interest
%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%
% G - Number for the governates of interest {'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
% Yt - Model output data for all governoates with non-zero incidence
% WID - incidence data for the governoates with non-zero incidence
% GNZI - the areas with non-zero incidence
% GTF - Gov in GNZI used in the fitting
%% Plot the figures

% Structure of the names of the gover
N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
NG=ceil(length(GNZI)./4); % length of the number of governates of interest
for ii=1:NG
    figure('units','normalized','outerposition',[0 0 1 1]); % opens a new figure
    for jj=1:4
        if(jj+4*(ii-1)<=length(GNZI))
            subplot('Position',[0.0708+rem(jj+1,2).*0.97/2,(0.9/2+0.163120567375887)-floor(jj/3).*0.9/2,0.897162184873949/2*0.95,0.793313069908819/2*0.82]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

            D=WID(GNZI(jj+4*(ii-1)),:); % Pull weekly incidence data for governonate

            MI=floor(log10(max(D))); % The max magnitude we want for the figure
            FF=ceil(max(D)./10^MI);
            NW=length(D); % Number of weeks
            scatter([1:NW],D,40,'k','filled'); hold on; % Plot data
            if(~isempty(find(GNZI(jj+4*(ii-1))==GNZI(GTF)))) % if included in the fitting
                plot([(1+(length(D)-length(Yt(jj+4*(ii-1),:)))):NW],Yt(jj+4*(ii-1),:),'k','LineWidth',2);  % Plot model predictions            
                plot([(1+(length(D)-length(Yt(jj+4*(ii-1),:)))):NW],Yt2(jj+4*(ii-1),:),'r','LineWidth',2); hold on; % Plot model predictions
            else
                plot([(1+(length(D)-length(Yt(jj+4*(ii-1),:)))):NW],Yt(jj+4*(ii-1),:),'-.','color',[0.5 0.5 0.5],'LineWidth',2);  % Plot model predictions            
                plot([(1+(length(D)-length(Yt(jj+4*(ii-1),:)))):NW],Yt2(jj+4*(ii-1),:),'-.','color',[0.7 0.2 0],'LineWidth',2); hold on; % Plot model predictions                
            end
            box off; % removes the outside box on the figure
            xlim([1 NW]); % sets the x-limits of our x -axis
            ylim([0 FF*10^MI]); %sets the y-limts of the y-axis
            % The size to separate the weeks in the x-label
            dW=7;
            % Set the X-tick labels to be Dates rather than numbers
            startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
            XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
            % changing the aspects of the axis for the the current figure 
            if(jj>2)
                set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:10^MI:FF*10^MI],'Fontsize',16,'XTickLabel',XTL);
            % Xlable of the figure
            xlabel('Week of report','Fontsize',24);
            else
                set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:10^MI:FF*10^MI],'Fontsize',16,'XTickLabel',{});
            end
            % Rotates the xticklabel 
            xtickangle(45) 
            % Sets the y-axis to not have 10^n
            ax=gca; % finds the current axis
            ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
            title(N(GNZI(jj+4*(ii-1))).G ,'Fontsize',28);
            if(rem(jj,2)==1)
            % ylable of the figure
            ylabel('Number of reported cases','Fontsize',24);
            end
        end
    end
end

end