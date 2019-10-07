
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% endPDS=0.8;
%% Names of the covariates
X=struct('N',{'Bias','Population Density','Health Facilities','WASH and Incidence','Incidence and Attacks','Incidence and Conflict','WASH, Rainfall, and Incidence','WASH and Rainfall','External Incidence','Attacks','Rebel Control','Conflict and Rainfall','Attack and Rainfall'});
%% Run the projection

PDS=0.8;

    atest=0.05;
    XUm=eye(13);

    NMR=length([1 2 3 4 9]);
    XUm=XUm([1 2 3 4 9],:);
    par=zeros(NMR,31);
    RSSv=zeros(NMR,1);
    CVE=zeros(NMR,1);
    for ii=1:NMR % Remove the covariate from the model
        [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,[-16.*ones(1,13) ones(1,7) 0 0 0 -16.*ones(1,8)]);
    end
    f=find(RSSv==min(RSSv));
    RSSv=RSSv(f);
    CVE=CVE(f);
    XUv=XUm(f,:);
    parv=par(f,:);
    [kv]=RetParameterPS(parv,XUv);

    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionNoCNoR(XUv,RSSv,kv,atest,PDS,par);

    save(['ForwardSelectionNoConflictNoRainfall-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);
    while(~isempty(XUr))
        par=[par; partemp];
        RSSv=[RSSv; RSSr];
        XUv=[XUv; XUr];
        parv=[parv; parr];
        kv=[kv; kr];
        CVE=[CVE; CVEr];
        [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionNoCNoR(XUv(end,:),RSSv(end),kv(end),atest,PDS,par);   
        save(['ForwardSelectionNoConflictNoRainfall-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);
    end
% delete pobj;