
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end

PDS=0.8;
%% Names of the covariates

% %% Run the projection
% NMR=13;
% XUm=eye(NMR);
% par=zeros(NMR,length([-16.*ones(1,NMR) ones(1,6) 0 0 0 -16.*ones(1,8)]));
% RSSv=zeros(NMR,1);
% CVE=zeros(NMR,1);
% for ii=1:NMR % Remove the covariate from the model
%     [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,[-16.*ones(1,NMR) ones(1,6) 0 0 0 -16.*ones(1,8)]);
% end
% f=find(CVE==min(CVE));
% RSSv=RSSv(f);
% CVE=CVE(f);
% XUv=XUm(f,:);
% parv=par(f,:);
% [kv]=RetParameterPS(parv,XUv);


atest=0;
load(['ForwardSelectionNoConflictNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);
atest=0;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,parv(end,:));

save(['ForwardSelection(SeedNoConflictNoRain)-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
while(~isempty(XUr))
    par=[partemp];
    RSSv=[RSSv; RSSr];
    CVE=[CVE; CVEr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,partemp);   
    save(['ForwardSelection(SeedNoConflictNoRain)-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
end
% delete pobj;