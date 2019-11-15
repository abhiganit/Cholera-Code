clear;
clc;

    XU=zeros(1,18);
    XU(1)=1;
    XU(2)=1;
%% Names of the covariates
X=struct('N',{'WASH & Food and Population Density','WASH & Food, Population Density, and Rebel Control','WASH & Food Security and Incidence','Population density, WASH & Food Security, and incidence','Health facilities, WASH & Food Security, and incidence','Rebel control','Targeted attacks and Incidence Indicator','Conflict and Incidence Indicator','Attack and Incidence Indicator','WASH, Incidence Indicator, and rainfall','WASH, rainfall, and incidence','Conflict, Incidence Indicator, and Rainfall','Targeted attack, Incidence Indicator and Rainfall','Attack, Incidence Indicator, and Rainfall','Wheat prices, Food security, and Incidence Indicator','Diesel prices, Food secuity, Incidence Indicator','Diesel prices, WASH, Incidence Indicator'});
lbps=[-32.*ones(1,length(XU)) zeros(1,length(XU)-7) -32.*ones(1,8) -32 -2 -32 -2 -32.*ones(1,7) -32 log10(0.9)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,length(XU)-7)  log10([ones(1,8) 20 1000 20 1000 120 120 120 120 120 1000 1 1000 exp(log(26/56)/(4*52))])]; % specify the upperbound for the parameters 

temppar=lbps+(ubps-lbps).*rand(size(lbps));


   %% Run the projection
   
   CF=[0;0];
      par=zeros(1,length(lbps));
        for jj=1:3
                [par,RSSv,CVE] =ProFittingGA(XU,0.8,CF,0,zeros(4,1),temppar);
        end