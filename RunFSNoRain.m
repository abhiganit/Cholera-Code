
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end

%% Names of the covariates

%% Run the projection
PDS=0.8;

%% Run the projection


atest=0;
load(['ForwardSelectionNoConflictNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);

atest=0;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionNoRain(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,parv(end,:));

save(['ForwardSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    CVE=[CVE; CVEr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelectionNoRain(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,partemp);   
    save(['ForwardSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(100*atest) '.mat']);
end
% delete pobj;