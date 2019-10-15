
if(isempty(gcp('nocreate')))
    pobj=parpool(20);
end
PDS=0.8;
%% Names of the covariates
X=struct('N',{'Bias','Population Density','Health Facilities','WASH and Incidence','Incidence and Attacks','Incidence and Conflict','WASH, Rainfall, and Incidence','Rainfall and Incidence','WASH, Pop. Desnity, and Incidence','External incidence IDP','Rebel Control','Conflict, Incidence, and Rainfall','Attack, Incidence and Rainfall'});
%% Run the projection
NMR=13;
XUm=eye(NMR);
par=zeros(NMR,length([-16.*ones(1,NMR) ones(1,6) 0 0 0 -16.*ones(1,8)]));
RSSv=zeros(NMR,1);
CVE=zeros(NMR,1);
for ii=1:NMR % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,[-16.*ones(1,NMR) ones(1,6) 0 0 0 -16.*ones(1,8)]);
end
f=find(RSSv==min(RSSv));
RSSv=RSSv(f);
CVE=CVE(f);
XUv=XUm(f,:);
parv=par(f,:);
[kv]=RetParameterPS(parv,XUv);

atest=0.05;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS,par);

save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
while(~isempty(XUr))
    par=[par;partemp];
    RSSv=[RSSv; RSSr];
    CVE=[CVE; CVEr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS,partemp);   
    save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
end
delete pobj;