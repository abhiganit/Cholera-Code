
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end
PDS=0.8;
%% Names of the covariates
X=struct('N',{'Bias','Population Density','Health Facilities','WASH and Incidence','Incidence and Attacks','Incidence and Conflict','Cumulative Attacks, Rainfall, and Incidence','Cumulative Attacks and Rainfall','External Incidence','Attacks','Rebel Control','Conflict and Rainfall','Attack and Rainfall'});
%% Run the projection

NMR=13;
XUm=eye(NMR);
par=zeros(NMR,31);
RSSv=zeros(NMR,1);
for ii=1:NMR % Remove the covariate from the model
    [par(ii,:),~,RSSv(ii)] = ProFittingMGA(XUm(ii,:),PDS,[],0,0,0,[-16.*ones(1,NMR) 10^(-3).*ones(1,10)  -16.*ones(1,8)]);
end
f=find(RSSv==min(RSSv));
RSSv=RSSv(f);
XUv=XUm(f,:);
parv=par(f,:);
[kv]=RetParameterPS(parv,XUv);

atest=0.05;
[XUr,RSSr,parr,kr] = ForwardSelectionMulti(XUv,RSSv,kv,atest,PDS,parv);

save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-MultiObj.mat']);
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = ForwardSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS,parv(end,:));   
    save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '-MultiObj.mat']);
end
% delete pobj;