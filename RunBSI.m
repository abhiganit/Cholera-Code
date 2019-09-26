if(isempty(gcp('nocreate')))
    pobj=parpool(20);
end
PDS=0.9;
XUv=[1 1 1 1 0 0 1 1 1 0 0];
[parv,~,RSSv]=ProFittingGA(XUv,PDS,[],0,0,0);
[kv]=RetParameterGA(parv,XUv);

atest=0.01;
[XUr,RSSr,parr,kr] = BackwardsSelection(XUv,RSSv,kv,atest,PDS);

save('BackwardsSelection.mat');
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = BackwardsSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS);   
    save(['BackwardsSelectionIncidence-PercentDataSet=' num2str(PDS*100) '.mat']);
end
delete pobj;