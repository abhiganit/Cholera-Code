clear;
clc;
%% Calibrates the saturation function only to imporove the fitting
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});
 
load('Combo.mat');
for nn=0:32      
    if(nn==0)
        load(['LSE-Fit-Vaccination-Null.mat'])
        RSSvt=RSSv;
        [par,fval_c]=ProFitting_sigma(XU,par);
        save(['LSE-Fit-Vaccination-Null-Calibrated_sigma.mat'],'par','XU','X','fval_c');          
    else        
        indx=INC{nn};
        load(['LSE-Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat'])
        RSSvt=RSSv;
        [par,fval_c]=ProFitting_sigma(XU,par);
        save(['LSE-Fit-Vaccination-IncidenceperCapita' C(indx).N '-Calibrated_sigma.mat'],'par','XU','X','fval_c');        
    end
end
   
 