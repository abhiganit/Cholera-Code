 %% Names of the covariates
X=struct('N',{'WASH & Food','Population Density','Health Facilities','WASH & Food Security and Incidence','Population density, WASH & Food Security, and incidence','Health facilities, WASH & Food Security, and incidence','Rebel control','Targeted attacks and Incidence Indicator','Conflict and Incidence Indicator','Attack and Incidence Indicator','WASH, Incidence Indicator, and rainfall','WASH, rainfall, and incidence','Conflict, Incidence Indicator, and Rainfall','Targeted attack, Incidence Indicator and Rainfall','Attack, Incidence Indicator, and Rainfall','Wheat prices, Food security, and Incidence Indicator','Diesel prices, Food secuity, Incidence Indicator','Diesel prices, WASH, Incidence Indicator'});
lbps=[-32.*ones(1,18) zeros(1,18-7) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,9) -32.*ones(1,3)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ubps=[ 5.*ones(1,18) ones(1,18-7) ones(1,7) log10([ones(1,8) 20 3 20 3 120 120 120 120 120 1000 1000 1000 1 1000 1000 1])];
temppar=lbps+(ubps-lbps).*rand(size(lbps));

for cvii=11:-1:1
    PDS=0.4+0.05.*(cvii-1);

   %% Run the projection
    NMR=length([1:6]);
    XUm=eye(18);
    XUm=XUm([1:6],:);
    
    par=zeros(NMR,length(lbps));
    RSSv=zeros(NMR,1);
    CVE=zeros(NMR,1);
    for ii=1:NMR % Remove the covariate from the model
        [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,temppar);
    end
    f=find(CVE==min(CVE));
    RSSv=RSSv(f);
    CVE=CVE(f);
    XUv=XUm(f,:);
    parv=par(f,:);
    [kv]=RetParameterPS(parv,XUv);
    IndxAllow=[1:6];
    gg=find(XUv(end,IndxAllow)==0);
    IndxAllow=IndxAllow(gg);
    atest=0;
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,par,IndxAllow);

    save(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
    while(~isempty(XUr))
        par=[partemp];
        RSSv=[RSSv; RSSr];
        CVE=[CVE; CVEr];
        XUv=[XUv; XUr];
        parv=[parv; parr];
        kv=[kv; kr];
        gg=find(XUv(end,IndxAllow)==0);
        IndxAllow=IndxAllow(gg);
        [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,partemp,IndxAllow);   
        save(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
    end
end