
% if(isempty(gcp('nocreate')))
%     pobj=parpool(20);
% end

% 
% %% Names of the covariates
% X=struct('N',{'Bias','Population Density','Health Facilities','WASH and Incidence','Population density, WASH, and incidence','Health facilities, WASH, and incidence','Rebel control','Targeted attacks and Incidence','Conflict and incidence','Attack and Incidence','WASH, Incidence, and rainfall','WASH and rainfall','Conflict, Incidence, and Rainfall','Targeted attack, Incidence and Rainfall','Attack, Incidence, and Rainfall'});
% %% Run the projection
% NMR=15;
% XUm=eye(15);
% par=zeros(NMR,47);
% RSSv=zeros(NMR,1);
% CVE=zeros(NMR,1);
% for ii=1:NMR % Remove the covariate from the model
%     [par(ii,:),RSSv(ii),CVE(ii)] = ProFittingGA(XUm(ii,:),PDS,[-32.*ones(1,15) zeros(1,8) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,5)]);
% end
% f=find(CVE==min(CVE));
% RSSv=RSSv(f);
% CVE=CVE(f);
% XUv=XUm(f,:);
% parv=par(f,:);
% [kv]=RetParameterPS(parv,XUv);
% 
PDS=0.8;
atest=0;

load(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
[XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,parv(end,:));

save(['ForwardSelection-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
while(~isempty(XUr))
    par=[partemp];
    RSSv=[RSSv; RSSr];
    CVE=[CVE; CVEr];
    XUv=[XUv; XUr];
    parv=[parv; parr];
    kv=[kv; kr];
    [XUr,RSSr,CVEr,parr,kr,partemp] = ForwardSelection(XUv(end,:),RSSv(end),CVE(end),kv(end),atest,PDS,partemp);   
    save(['ForwardSelection-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
end
% delete pobj;