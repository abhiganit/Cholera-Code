function GovFigureIncidence(G,Yt,WID,GNZI)
% Plots the incidence data and model prediction for the area of interest
%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%
% G - Number for the governates of interest {'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
% Yt - Model output data for all governoates with non-zero incidence
% WID - incidence data for the governoates with non-zero incidence
% GNZI - the areas with non-zero incidence

%% Plot the figures

% Structure of the names of the gover
N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
NG=length(G); % length of the number of governates of interest
for ii=1:NG
    ff=find(G(ii)==GNZI); % see if the area is one with non-zero incidcen
    if(~isempty(ff))
        figure('units','normalized','outerposition',[0 0 1 1]); % opens a new figure
        subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

        D=WID(ff,:); % Pull weekly incidence data for governonate

        MI=ceil(log10(max(D))); % The max magnitude we want for the figure
        NW=length(D); % Number of weeks
        M=Yt(ff,:); % Pull weekly incidence from model for governonate
        scatter([1:NW],D,40,'k','filled'); hold on; % Plot data
        plot([(1+(length(D)-length(M))):NW],M,'b','LineWidth',2); hold off; % Plot model predictions

        box off; % removes the outside box on the figure
        xlim([1 NW]); % sets the x-limits of our x -axis
        ylim([0 10.^MI]); %sets the y-limts of the y-axis
        % The size to separate the weeks in the x-label
        dW=4;
        % Set the X-tick labels to be Dates rather than numbers
        startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
        XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
        % changing the aspects of the axis for the the current figure 
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:10^(MI-1):10.^MI],'Fontsize',16,'XTickLabel',XTL);
        % Rotates the xticklabel 
        xtickangle(45) 
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        title(N(G(ii)).G ,'Fontsize',28);
        % Xlable of the figure
        xlabel('Week of report','Fontsize',24);

        % ylable of the figure
        ylabel('Number of reported cases','Fontsize',24);
    else % area has zero incidence
        fprintf([ N(G(ii)).G ' does not have any reported incidence and was removed from the fitting process \n']); % print area with zero incidence
    end
end

end