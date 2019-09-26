% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end
PDS=1;
XUv=ones(1,11);
[parv,~,RSSv]=ProFittingGA(XUv,PDS,[],0,0,0,[-16.*ones(1,11) ones(1,5) 0 0 0 -16.*ones(1,8)]);
[kv,~]=RetParameterPS(parv,XUv);

atest=0.01;
[XUr,RSSr,parr,kr] = BackwardsSelection(XUv,RSSv,kv,atest,PDS);

save(['BackwardsSelection-PercentDataSet=' num2str(PDS*100) '.mat']);
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = BackwardsSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS);   
    save(['BackwardsSelection-PercentDataSet=' num2str(PDS*100) '.mat']);
end
delete pobj;