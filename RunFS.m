
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end

PDS=1;
%% Names of the covariates
Gov=9;
%% Names of the covariates
X=struct('N',{'Bias','Incidence','Incidence and taegeted attacks','Incidence and conflict','Incidence and attacks','Rainfall','Rainfall and conflict','Rainfall and targeted attacks','Rainfall and attacks'});
%% Run the projection
NMR=9;
XUm=eye(9);
par=zeros(NMR,length([-32.*ones(1,9) 10^(-16).*ones(1,8) 10^(-16).*ones(1,6) -32.*ones(1,6) -32.*ones(1,4) -32.*ones(1,4)]));
RSSv=zeros(NMR,1);
CVE=zeros(NMR,1);
for ii=1:NMR % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,Gov,[-32.*ones(1,9) 10^(-16).*ones(1,8) 10^(-16).*ones(1,6) -32.*ones(1,6) -32.*ones(1,4) -32.*ones(1,4)]);
end
f=find(RSSv==min(RSSv));
RSSv=RSSv(f);
CVE=CVE(f);
XUv=XUm(f,:);
parv=par(f,:);
[kv]=RetParameterPS(parv,XUv);

atest=0.05;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,Gov,parv(end,:));

save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-Gov=' num2str(Gov) '.mat']);
while(~isempty(XUr))
    par=[partemp];
    RSSv=[RSSv; RSSr];
    CVE=[CVE; CVEr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,Gov,partemp);   
    save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-Gov=' num2str(Gov) '.mat']);
end
% delete pobj;