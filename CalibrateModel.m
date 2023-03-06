clear;
clc;
%% Calibrates the saturation function only to imporove the fitting
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});
 
load('Combo.mat');
for nn=1:32             
    indx=INC{nn};
    load(['Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat'])
    RSSvt=RSSv;
    [par,DAR,RSSv]=ProFitting(XU,par);
    save(['Fit-Vaccination-IncidenceperCapita' C(indx).N '-CalibratedDAR.mat'],'par','RSSv','XU','X','DAR');        
end
   
 