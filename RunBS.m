if(isempty(gcp('nocreate')))
    pobj=parpool(20);
end
XUv=ones(1,11);
[parv,RSSv]=ProFittingGA(XUv,[],0,0,0);
[kv]=RetParameterGA(parv,XUv);

atest=0.01;
[XUr,RSSr,parr,kr] = BackwardsSelection(XUv,RSSv,kv,atest);

save('BackwardsSelection.mat');
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = BackwardsSelection(XUv(end,:),RSSv(end),kv(end),atest);   
    save('BackwardsSelection.mat');
end
delete pobj;