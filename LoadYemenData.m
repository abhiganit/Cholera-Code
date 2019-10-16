function [WI,Ctv,tA,Rtv,Mt,maxtau] = LoadYemenData
% Loads the Data needed to rum the regression model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WI - Weekly incidence (22x153)
% Ctv - Weekly level of conflict (22x153)
% tA - Weekly number of attacks (22x153)
% Rtv - Weekly level of rainfall (22x153)
% Mt - Water course connection (22x22)
% P - log of population density (22x1)
% RC - Rebel control (22x1)
% H - Health sites per 10000 (22x1)
% WPIN - transformation of density of WASH people in need (22x1)
% Et - External incidence due to IDP 
% GNZI - Gov with non-zero incidence
% maxtau - the maximum lag considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('Yemen_Gov_Incidence.mat'); % Incidence data
    S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

    % Record the names for the IDP calculation
    Sm=cell(length(S),1); % allocate space
    for ii=1:length(Sm)
    Sm{ii}=S(ii).ADM1_EN; % record name
    end

    WI=IData'; % Transpose the data set such that the number of areas is the row
    load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
    Ctv=GLevelConflict(ProC,S,153);
    load('Yemen_Air_Shelling.mat');
    Mt=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
    % Use attacks with "water in description"
    load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
    tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
    load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
    Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified
    % W = shaperead([ pwd '\ShapeFile\Wadies.shp']); % Shape file for Yemen water course
    % Mt=WaterCourseConnection(W,S); % Calculate the contact matrix for the water course
    % Load population density
    maxtau=4; % The maximum lag allowed for the model


end

