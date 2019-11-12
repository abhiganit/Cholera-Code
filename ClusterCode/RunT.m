   pobj=parpool(20);
%% Names of the covariates
X=struct('N',{'WASH & Food','Population Density','Health Facilities','WASH & Food Security and Incidence','Population density, WASH & Food Security, and incidence','Health facilities, WASH & Food Security, and incidence','Rebel control','Targeted attacks and Incidence Indicator','Conflict and Incidence Indicator','Attack and Incidence Indicator','WASH, Incidence Indicator, and rainfall','WASH, rainfall, and incidence','Conflict, Incidence Indicator, and Rainfall','Targeted attack, Incidence Indicator and Rainfall','Attack, Incidence Indicator, and Rainfall','Wheat prices, Food security, and Incidence Indicator','Diesel prices, Food secuity, Incidence Indicator','Diesel prices, WASH, Incidence Indicator'});
lbps=[-32.*ones(1,18) zeros(1,18-7) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,9) -32.*ones(1,3)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,18) ones(1,18-7) ones(1,7) log10([ones(1,8) 20 1000 20 1000 120 120 120 120 120 1000 1000 1000 1 1000 1000 1])];
temppar=lbps+(ubps-lbps).*rand(size(lbps));
XU=zeros(1,18);
XU(1)=1;
XU(4)=1;
XU(7)=1;
XU(8)=1;
XU(9)=1;
XU(10)=1;

XU(16)=1;
XU(17)=1;
XU(18)=1;

par=zeros(11,length(lbps));
RSSv=zeros(11,1);
CVE=zeros(11,1);
PDS=0.4+0.05.*[0:10];
parfor cvii=1:11    

   %% Run the projection
        [par(cvii,:),RSSv(cvii),CVE(cvii)] = ProFittingGA(XU,PDS(cvii),temppar);

end

    save(['Fit-Vaccination-PercentDataVARY.mat']);
    
    
delete pobj;