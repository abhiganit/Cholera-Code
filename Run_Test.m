clear;
clc;

    XU=zeros(1,18);
XU(4)=1;
XU(7)=1;
XU(8)=1;
XU(9)=1;
XU(10)=1;

XU(16)=1;
XU(17)=1;
XU(18)=1;
%% Names of the covariates
X=struct('N',{'WASH & Food','Population Density','Health Facilities','WASH & Food Security and Incidence','Population density, WASH & Food Security, and incidence','Health facilities, WASH & Food Security, and incidence','Rebel control','Targeted attacks and Incidence Indicator','Conflict and Incidence Indicator','Attack and Incidence Indicator','WASH, Incidence Indicator, and rainfall','WASH, rainfall, and incidence','Conflict, Incidence Indicator, and Rainfall','Targeted attack, Incidence Indicator and Rainfall','Attack, Incidence Indicator, and Rainfall','Wheat prices, Food security, and Incidence Indicator','Diesel prices, Food secuity, Incidence Indicator','Diesel prices, WASH, Incidence Indicator'});
lbps=[-32.*ones(1,length(XU)) zeros(1,length(XU)-7) -32.*ones(1,8) -32 -2 -32 -2 -32.*ones(1,7) -32 log10(0.9)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,length(XU)) ones(1,length(XU)-7)  log10([ones(1,8) 20 1000 20 1000 120 120 120 120 120 1000 1 1000 exp(log(26/56)/(4*52))])]; % specify the upperbound for the parameters 

temppar=lbps+(ubps-lbps).*rand(size(lbps));


   %% Run the projection
   
   CF=[0 1 2;0 1 2];
      par=zeros(3,length(lbps));
      RSSv=zeros(3,1);
      CVE=zeros(3,1);
        for jj=1:3
           parfor ii=1:3
                [par(ii,:),RSSv(ii),CVE(ii)] =ProFittingGA(XU,0.8,CF(:,ii),0,zeros(4,1),temppar);
           end

            temppar=par;
            for kk=1:3
                for ii=1:3               
                    [lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par(ii,:));
                    % Take paramters ii and plug into function CF for kk
                    Rtemp=OFuncProPS(pars,CF(:,kk),0,zeros(4,1),WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW));
                    if(Rtemp<RSSv(kk))
                       temppar(kk,:)=par(ii,:); 
                       RSSv(kk)=RSSv(ii);
                       CVE(kk)=CVE(ii);
                    end
                end
            end
            par=temppar;
            temppar=unique(temppar,'rows');
            save(['Fit-Vaccination-PercentData=80.mat'],'par','RSSv','CVE','XU','X','CF');
        end
