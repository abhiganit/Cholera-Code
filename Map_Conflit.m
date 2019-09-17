%% Plots the heat map of the conflict and produces cross-hatch in seprate figure for specified areas

%% Construction of the heat map
close all; % Close all figures
load('Conflict_Yemen_Time_Location_Project_Forward'); % load the conflict data
Loc=ProC(:,2:3); % remove the time component from the conflic tdata
 h=histogram2(Loc(:,2),Loc(:,1),[75 51]); % construct a 2d histogram of where the events took place
 X=(h.XBinEdges(2:end)+h.XBinEdges(1:end-1))./2; % determine the center of the bins in x irection
 Y=(h.YBinEdges(2:end)+h.YBinEdges(1:end-1))./2; % determine the center of hte bins in the y direction
 Z=log(h.Values+1)'; % condcut a log -transform on the evetns
 close all; % close historgam
 Xq=linspace(min(X),max(X),301); % refine the x axis
 Yq=linspace(min(Y),max(Y),201); % refine the y-axis
 [Xm,Ym]=meshgrid(Xq,Yq); % greates a matrix of the (x,y) values
 ZZ=interp2(X,Y,Z,Xm,Ym); % Interpolates the data based on the values from the histogram
 

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
figure('units','normalized','outerposition',[0 0 1 1]); % Crat a new figure
 contourf(Xq,Yq,ZZ,'LineStyle','none'); hold on; % plots a contour map with no line
 load('HeatMap.mat','cmap'); % loads a custom color map
 colormap(cmap); % changes the color map of the figure
 mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on; % plot the map of yemen 
 

set(gca,'visible','off') % remove the asix from the figures
%print(gcf,'Conflict_Map_Yemen.png','-dpng','-r600'); % saves the conflict map

%% Produces the map with a cross-hatch in the specified areas
figure('units','normalized','outerposition',[0 0 1 1]); % opens a new figure

SC = shaperead([ pwd '\ShapeFile\yem_admbnda_adm0_govyem_mola_20181102.shp']); % the shape file of only yemen

 mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on; % plots the lines of the governorates oand no face colour
CrossHatchMap(min(SC.X),max(SC.X),min(SC.Y),max(SC.Y),0.005,[0 0 0],SC); hold on % produces a cross hatch of yemen as a whole
% Define a structure of the names for the different Governorate
G=struct('Name',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
% Determines the number of Governorate from the structure above
NG=length(G); % the number of governorates 
% A vector specifying the the area of Yemen the Governorate belongs to
% (Index starts at 1 as we will be using these numbers to call the colours
% we nedd when creating the figure

AG=[3;3;1;1;1;2;4;2;2;2;1;4;2;1;3;2;1;2;2;4;4;1]; % the segragation of th different governorates

Sm=cell(length(S),1); % cell array for the names of the areas
for ii=1:length(Sm)
    Sm{ii}=S(ii).ADM1_EN; % record the names specified in the shape file
end
for ii=1:NG % go through the governorates
    ff=find(contains(Sm,G(ii).Name)); % see if the name in shape file matchs that in structures
    if(AG(ii)<=2) % if the area is specifeid as 1 or two then plot with white face color
        mapshow(S(ff),'FaceColor','w','Edgecolor',[0 0 0],'LineWidth',1);hold on %plot with white face color
    end
end
set(gca,'visible','off') % removes the axis
% print(gcf,'Conflict_Map_Yemen_Primary_Outbreak.png','-dpng','-r600'); % saves the figure