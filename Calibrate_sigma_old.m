clear;
clc;

%% Calibrates the saturation function only to imporove the fitting
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});
 
load('Combo.mat');
for nn=1:32      
        indx=INC{nn};
        
        XU=zeros(1,32);
        XU([21:32])=1; % Turn on rainfall, temprature and incidence
        for ii=1:length(indx)
            XU([1:4]+4.*(indx(ii)-1))=1;     
        end
        load(['Fit-Vaccination-IncidenceperCapita' C(indx).N '-Rain.mat'],'par')
        par_t=[par(1:24) -32.*ones(1,4) par(25:28) par(29:36) -32 par(37:39) -3];
        [par,fval_c]=ProFitting_sigma(XU,par_t);
        save(['OLD-Fit-Vaccination-IncidenceperCapita' C(indx).N '-Rain-Calibrated_sigma.mat'],'par','XU','fval_c');        
end


