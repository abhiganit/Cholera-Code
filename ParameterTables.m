load('Combo.mat');
beta=zeros(64,28);
DA=zeros(64,1);
K=zeros(2,64);
n=zeros(2,64);
dV=zeros(2,64);
KP=zeros(2,64);
r0=zeros(64,1);
KV=zeros(64,1);
w=zeros(64,1);
DARv=zeros(64,1);
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
for ii=1:64
    load(['Fit-Vaccination-IncidenceperCapita' C(INC{ii}).N '-CalibratedDAR.mat']);
    [~,beta(ii,:),tau,DB,DA(ii),K(:,ii),n(:,ii),KP(:,ii),KV(ii),dV(:,ii),r0(ii),~,w(ii)]=RetParameterPS(par,XU,CF,4);
    DARv(ii)=DAR;
end
K=K';
n=n';
KP=KP';
dV=dV';
T=table(beta,DA,K,n,KP,KV,dV,r0,DARv,w);
 
 writetable(T,'ModelParameter.csv','Delimiter',',');