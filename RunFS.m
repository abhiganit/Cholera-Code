
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end
PDS=1;

%% Run the projection

NMR=11;
XUm=eye(NMR);
par=zeros(NMR,27);
RSSv=zeros(NMR,1);
for ii=1:NMR % Remove the covariate from the model
    [par(ii,:),~,RSSv(ii)] = ProFittingGA(XUm(ii,:),PDS,[],0,0,0,[-16.*ones(1,11) ones(1,5) 0 0 0 -16.*ones(1,8)]);
end
f=find(RSSv==min(RSSv));
RSSv=RSSv(f);
XUv=XUm(f,:);
parv=par(f,:);
[kv]=RetParameterPS(parv,XUv);

atest=0.01;
[XUr,RSSr,parr,kr] = ForwardSelection(XUv,RSSv,kv,atest,PDS,parv);

save('BackwardsSelection.mat');
while(~isempty(XUr))
    RSSv=[RSSv; RSSr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,parr,kr] = ForwardSelection(XUv(end,:),RSSv(end),kv(end),atest,PDS,parv(end,:));   
    save(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '.mat']);
end
% delete pobj;