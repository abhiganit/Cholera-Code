function Ct= GLevelConflict(C,S,NW)
% Ct= GLevelConflict(C,S,NW) constructs a matrix of the of conflict/attacks
% agonf the areas over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C - a matrix of size Nx3 where N is the numebr of event, the first column
% is the week (i.e. 1 indicates the first week (i.e. week of Oct 3 2016),
% the second column is the latitude and the third column is the longitude
% S - the shape file for the area of interest (M areas in the shape file)
% NW - the number of weeks to be included for the matrix Ct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ct - the number of events that happened in an area over various weeks
% (size MxNW )

%% Construct the matrix
M=length(S); % the number of areas in the shape file
f=find(C(:,1)<=NW); % Only want to include the conflict that is within the time span we want
C=C(f,:); % truncate to the weeks that we want for there is no index error for events outside the time of interest
Ct=zeros(M,NW); % Construct the matrix for the number of events
for rr=1:M % go throough all M areas in the shape file
    in=inpolygon(C(:,3),C(:,2),S(rr).X,S(rr).Y);  % Returns whether the points in C are in the polygon spcified by S 
    f=find(in>0); % Find if points are in territory
    for jj=1:length(f)
       Ct(rr,C(f(jj),1))= Ct(rr,C(f(jj),1))+1; % Adds the number of events in area rr to the specified week of the event C(f(jj),1) for the points
    end
end

end

