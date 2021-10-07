clear;
clc;
%% Calibrates the saturation function only to imporove the fitting
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
 
load('Combo.mat');
for nn=1:64             
    indx=INC{nn};
    load(['Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat'])
    RSSvt=RSSv;
    [par,DAR,RSSv]=ProFitting(XU,CF,RF,par);
    save(['Fit-Vaccination-IncidenceperCapita' C(indx).N '-CalibratedDAR.mat'],'par','RSSv','XU','CF','RF','X','DAR');        
end
   
 