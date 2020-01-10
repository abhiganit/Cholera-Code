
load('ForwardSelectionNoRainNoConflict-alpha=0-PercentData=65.mat');
IndxAllow=[1:10 16 17 18];
gg=find(XUv(end,IndxAllow)==0);
IndxAllow=IndxAllow(gg);
atest=0;
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,par,IndxAllow);

save(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
while(~isempty(XUr))
    par=[partemp];
    RSSv=[RSSv; RSSr];
    CVE=[CVE; CVEr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    gg=find(XUv(end,IndxAllow)==0);
    IndxAllow=IndxAllow(gg);
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,partemp,IndxAllow);   
    save(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
end
