if(isempty(gcp('nocreate')))
    pobj=parpool(20);
end
XUv=[1 1 1 1 0 0 1 1 1 0 0];
[parv,RSSv]=ProFittingGA(XUv,[],0,0,0);
[kv]=RetParameterGA(parv,XUv);

atest=0.01;
[XUr,RSSr,parr,kr] = BackwardsSelection(XUv,RSSv,kv,atest);
  
save('BackwardsSelection_Incidence.mat');
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = BackwardsSelection(XUv(end,:),RSSv(end),kv(end),atest);   
    save('BackwardsSelection_Incidence.mat');
end
delete pobj;