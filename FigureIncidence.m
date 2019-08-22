%% This script creates four panels to show the weekly incidence among the different Governorate and different areas in Yemen
% Close all previous windows
close all;
% clear all previous data
clear; 
% Define a structure of the names for the different Governorate
G=struct('Name',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
% Determines the number of Governorate from the structure above
NG=length(G);
% A vector specifying the the area of Yemen the Governorate belongs to
% (Index starts at 1 as we will be using these numbers to call the colours
% we nedd when creating the figure
AG=[3;3;1;1;1;2;4;2;2;2;1;4;2;1;3;2;1;2;2;4;4;1];

% Specify the colours/shades for the Governorate
% Colour range for area 1
CM=hex2rgb(['#A11F0C';'#FFA577']); % This functino converts hexdecimial colors to and rgb vecotr: The colorurs here are red
% Colour range for area 2
CM=[CM;hex2rgb(['#C99E10';'#E6D72A'])];  % This functino converts hexdecimial colors to and rgb vecotr: The colorurs here are yellow
% Colour range for area 3
CM=[CM;hex2rgb(['#515F37';'#c7e9c0'])];  % This functino converts hexdecimial colors to and rgb vecotr: The colorurs here are green
% Colour range for area 4
CM=[CM;hex2rgb(['#3f007d';'#dadaeb'])];  % This functino converts hexdecimial colors to and rgb vecotr: The colorurs here are purple/brown
% Add other colours if we need more areas defined

% Colour matrix to be used in the construction of the figure
C=zeros(NG,3); % We have NG Governorates and the rgb vector is 3 numbers
% Find how many Governorate in area ii
for ii=1:max(AG)
    h=length(find(AG==ii)); % Find and count the Governorates in area ii 
    if(h>0) % only do the next steps if there is at least one Governorates in area ii 
        idx = linspace(0,1,h);    % split the interval from 0 to 1 into a vector of size h that is evenly 
        temp = interp1([0 1],CM([1:2]+2.*(ii-1),:),idx); % Steps throught the matrix 2 rows at a time that we created for the different areas in Yemen
        CD=length(find(AG<ii)); % This find the number we have already completed
        for jj=1:h % Adding the number of colours to the C matrix
           C(CD+jj,:)= temp(jj,:); % puts the jj index from temp into the proper index in the colour matrix in position f(jj)
        end
    end
end

OG=zeros(NG,1); % Specifys the odering for the different governanents based on the area such that the reds are with reds, etc..
for ii=1:max(AG)
  f=find(AG==ii); % find which area it is in
  g=find(OG==0,1); % finds the first index in OG that has not been filled yet
  OG([0:length(f)-1]+g)=f; % finds the index of choice to re-arrange the incidence data such that our colours mathc up correctly
end

load('Yemen_Gov_Incidence.mat'); % load governorate incidence data

% Creates a figure that we can later refer to as f1
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
NW=length(IData(:,1)); % Determine the number of weeks that we have for the Governorates IData is NWxNG in size
hb=area([1:NW],IData(:,OG),'LineStyle','none'); % this plots a stacked area graph based on the weekly incidence data from IData. We remove the line from the grapgh with 'LineStyle','none'
% This code here changes the colours of the area graph to the ones that we
% want
for mm=1:NG
    hb(mm).FaceColor=C(mm,:); % Changing the colour for each of the specified Governorates
end
box off; % removes the outside box on the figure
xlim([1 NW]); % sets the x-limits of our x -axis
ylim([0 55000]); %sets the y-limts of the y-axis
% The size to separate the weeks in the x-label
dW=4;
% Set the X-tick labels to be Dates rather than numbers
startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:5000:55000],'Fontsize',16,'XTickLabel',XTL);
% Rotates the xticklabel 
xtickangle(45) 
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n

% Xlable of the figure
xlabel('Week of report','Fontsize',24);

% ylable of the figure
ylabel('Number of reported cases','Fontsize',24);

% Puts the legend of the Governorates through through coloured text
for ii=1:NG
   text(1.25,max(ylim)-1750*(ii-1)-500,G(OG(ii)).Name,'color',C(ii,:),'Fontsize',16);
end

% Creates a figure that we can later refer to as f1
f1=figure('units','normalized','outerposition',[0 0 1 1]);

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % load the shape file
Sm=cell(length(S),1); % cell array for the names of areas in shape file
for ii=1:length(Sm)
    Sm{ii}=S(ii).ADM1_EN; % record the names of the shape file
end

mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on % plots the lines of the borders
for ii=1:NG % shade the area of the govneroate
    ff=find(contains(Sm,G(OG(ii)).Name)); % see if the name is in the shape file based on the ordering OG
    mapshow(S(ff),'FaceColor',C(ii,:),'Edgecolor',[0 0 0],'LineWidth',1);hold on % shape the area in with color C(ii,:)
end
xl=xlim; % saves the xlim to xl
yl=ylim;% saves the ylim to yl
% Produces a cross hatch for polygon 12
CrossHatchMap(min(S(12).X),max(S(12).X),min(S(12).Y),max(S(12).Y),0.01,[1 1 1],S(12))
% Produces a cross hatch for polygon 21
CrossHatchMap(min(S(21).X),max(S(21).X),min(S(21).Y),max(S(21).Y),0.025,[0 0 0],S(21))
xlim(xl); % set xlimits to xl
ylim(yl); % set ylimits to yl
set(gca,'visible','off') % remove the axis from the map 