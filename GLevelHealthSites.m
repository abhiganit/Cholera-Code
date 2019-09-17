function Ht= GLevelHealthSites(H,S)
%GLEVELHEALTHSITES reutrns the number of health sites in each polygon ofthe
%shape file S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H - Shapefile for the health sites
% S - Shapefile for the country
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ht- the vector of the number of health sites

%% Construct the matrix
M=length(S); % the number of areas in the shape file
Ht=zeros(M,1); % Construct the matrix for the number of events
HX=zeros(length(H),1);
HY=zeros(length(H),1);
for ii=1:length(H)
    HX(ii)=H(ii).X;
    HY(ii)=H(ii).Y;
end
for rr=1:length(Ht)
    in=inpolygon(HX,HY,S(rr).X,S(rr).Y);  % Returns whether the points in C are in the polygon spcified by S 
    f=find(in>0); % Find if points are in territory
    Ht(rr)=length(f);
end
end

