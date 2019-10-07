
if(isempty(gcp('nocreate')))
    pobj=parpool(20);
end
% Percentage of the data set wanting to fit
PDS=0.8;
%% Names of the covariates
X=struct('N',{'Bias','Population Density','Health Facilities','WASH and Incidence','Incidence and Attacks','Incidence and Conflict','Rainfall, and Incidence','Rainfall','External Incidence','Attacks','Rebel Control','Conflict and Rainfall','Attack and Rainfall'});
%% Run the projection

NMR=13; % Number of models to run
XUm=eye(NMR); % construct the covairates to include
par=zeros(NMR,31); % aloocate memory par
RSSv=zeros(NMR,1); % aloocate memory par
CVE=zeros(NMR,1);
for ii=1:NMR % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingMGA(XUm(ii,:),PDS,[-16.*ones(1,NMR) 10^(-3).*ones(1,10)  -16.*ones(1,8)]);
end
f=find(RSSv==min(RSSv));
RSSv=RSSv(f);
XUv=XUm(f,:);
parv=par(f,:);
CVE=CVE(f);
[kv]=RetParameterPS(parv,XUv);

atest=0.01;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionMulti(XUv,RSSv,kv,atest,PDS,par);

save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-MultiObj.mat']);
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    CVE=[CVE;CVEr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionMulti(XUv(end,:),RSSv(end),kv(end),atest,PDS,partemp);   
    save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-MultiObj.mat']);
end
delete pobj;