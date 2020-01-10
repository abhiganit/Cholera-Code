clear;
clc;
%% Calibrates the saturation function only to imporove the fitting
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
 

indx=[1 2 3 5];
    for yy=1:4 
        for ss=yy:4
            load(['Fit-Vaccination-PercentData=80-IncidenceperCapita-Diesel-Rain' C(unique([indx(ss) indx(yy)])).N '.mat']);
            [par,DAR,RSSv,CVE]=ProFitting(XU,0.8,CF,RF,par);
            save(['Fit-Vaccination-PercentData=80-IncidenceperCapita-Diesel-Rain' C(unique([indx(ss) indx(yy)])).N '-Calibrate-DAR.mat'],'par','RSSv','CVE','XU','CF','RF','X','DAR');
        end
    end
   
 