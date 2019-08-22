function CrossHatchMap(minx,maxx,miny,maxy,ds,c,S)
%% CrossHatchMap(minx,maxx,miny,maxy,ds)
% Constructs a cross hatch for the underlay of the map to plot the impact
% of attacks on Reff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minx - the minimum x value of the map
% maxx - the mximum x value of the map
% miny - the minimum y value of the map
% maxy - the mximum y value of the map
% ds - the spacing between the lines
% c - color of the lines
% S - the area you want the cross hatch contained in

xx=linspace(minx+miny-maxy,maxx,round(1/ds)); % Specify the points where we want line

CH1=cell(length(xx),1);
CH2=cell(length(xx),1);
SS=501; % the refinement of the areas
dx=linspace(0,maxy-miny,SS); % the dsiatance between the points 
dx=dx(2)-dx(1); % the unifrom distance

%% Determine the points in the polygon
parfor ii=1:length(xx) % loop in parallel through the values of xx
    xp=linspace(xx(ii),maxy-miny+xx(ii),SS); % the x-values for the point starting at xx(ii)
    yp=(linspace(xx(ii),maxy-miny+xx(ii),SS)-xx(ii))+miny; % the y-values for the point starting at xx(ii) for the postive slope
    f=inpolygon(xp,yp,S.X,S.Y); % see which points lie in the specified shape file
    CH1{ii}=[xp(f); yp(f)]; % record the points that lie in the polygon

   yp=-(linspace(xx(ii),maxy-miny+xx(ii),SS)-xx(ii))+maxy;   % the y-values for the point starting at xx(ii) for the negative slope
   f=inpolygon(xp,yp,S.X,S.Y); % see which points lie in the specified shape file

    CH2{ii}=[xp(f); yp(f)]; % record the points that lie in the polygon
end

%% Plot the lines in the shape file
for ii=1:length(xx)   % go through all intersectino points
    P=CH1{ii}; % Points of interst to be ploted
    dP=P(1,2:end)-P(1,1:end-1); % see if there are any discontiuities in the points
    f=find(abs(dP-dx)>dx/2); % ensure that they are not larger than univorm being specified
    if(~isempty(f)) % if f is not empty then polot the lines in f
        if(f(1)>1) % Go the first point where there is the disconituty
            plot(P(1,1:(f(1)-1)),P(2,1:(f(1)-1)),'color',c,'LineWidth',0.5); hold on;  % plot the line
            lp=length(P(1,f(1):end)); % adjust the length of the remaining points to be plotted
            P=P(:,f(1):end); % Truncate P from the point where discontunity starts
        else
            lp=length(P(1,(f(1)+1):end)); % adjust the length of the remaining points to be plotted, we skip as there is a disconuity
            P=P(:,(f(1)+1):end); % Truncate P from the point where discontunity starts
        end
    else   % if empty then plot all points         
        plot(P(1,1:end),P(2,1:end),'color',c,'LineWidth',0.5); hold on; 
        lp=0; % no points remaining
    end
   while(lp>1) % contiunue until on one point remianing
        dP=P(1,2:end)-P(1,1:end-1); % see if there are any discontiuities in the points
        f=find(abs(dP-dx)>dx/2); % ensure that they are not larger than univorm being specified
        if(~isempty(f)) % if f is not empty then polot the lines in f
            if(f(1)>1) % Go the first point where there is the disconituty
                plot(P(1,1:(f(1)-1)),P(2,1:(f(1)-1)),'color',c,'LineWidth',0.5); hold on;  % plot the line
                lp=length(P(1,f(1):end)); % adjust the length of the remaining points to be plotted
                P=P(:,f(1):end); % Truncate P from the point where discontunity starts
            else
                lp=length(P(1,(f(1)+1):end)); % adjust the length of the remaining points to be plotted, we skip as there is a disconuity
                P=P(:,(f(1)+1):end); % Truncate P from the point where discontunity starts
            end
        else   % if empty then plot all points         
            plot(P(1,1:end),P(2,1:end),'color',c,'LineWidth',0.5); hold on; 
            lp=0; % no points remaining
        end
   end
    P=CH2{ii}; % Points to be plotted
    dP=P(1,2:end)-P(1,1:end-1); % see if there are any discontiuities in the points
    f=find(abs(dP-dx)>dx/2); % ensure that they are not larger than univorm being specified
    if(~isempty(f)) % if f is not empty then polot the lines in f
        if(f(1)>1) % Go the first point where there is the disconituty
            plot(P(1,1:(f(1)-1)),P(2,1:(f(1)-1)),'color',c,'LineWidth',0.5); hold on;  % plot the line
            lp=length(P(1,f(1):end)); % adjust the length of the remaining points to be plotted
            P=P(:,f(1):end); % Truncate P from the point where discontunity starts
        else
            lp=length(P(1,(f(1)+1):end)); % adjust the length of the remaining points to be plotted, we skip as there is a disconuity
            P=P(:,(f(1)+1):end); % Truncate P from the point where discontunity starts
        end
    else   % if empty then plot all points         
        plot(P(1,1:end),P(2,1:end),'color',c,'LineWidth',0.5); hold on; 
        lp=0; % no points remaining
    end
   while(lp>1)
        dP=P(1,2:end)-P(1,1:end-1); % see if there are any discontiuities in the points
        f=find(abs(dP-dx)>dx/2); % ensure that they are not larger than univorm being specified
        if(~isempty(f)) % if f is not empty then polot the lines in f
            if(f(1)>1) % Go the first point where there is the disconituty
                plot(P(1,1:(f(1)-1)),P(2,1:(f(1)-1)),'color',c,'LineWidth',0.5); hold on;  % plot the line
                lp=length(P(1,f(1):end)); % adjust the length of the remaining points to be plotted
                P=P(:,f(1):end); % Truncate P from the point where discontunity starts
            else
                lp=length(P(1,(f(1)+1):end)); % adjust the length of the remaining points to be plotted, we skip as there is a disconuity
                P=P(:,(f(1)+1):end); % Truncate P from the point where discontunity starts
            end
        else   % if empty then plot all points         
            plot(P(1,1:end),P(2,1:end),'color',c,'LineWidth',0.5); hold on; 
            lp=0; % no points remaining
        end
   end
end
xlim([minx maxx]); % set the limits for the x-axis
ylim([miny maxy]); % set the limits for the y-axis
end

