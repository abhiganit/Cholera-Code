% Evaluates the parameters based on the leave 3 out cross-validation


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

C = combnk([1:21],3);

NN=length(C(:,1));
IG=zeros(NN,length(GNZI));
TotA=sum(tA,2);
TotA=find(TotA>0);
GV=find(GV==1);
for kk=1:NN
    % Gov included in the fitting
    GTF=setdiff([1:21],C(kk,:));
    IG(kk,GNZI(GTF))=1;
end

clearvars -except IG


load('KFold-Vaccination-PercentData=80-IncidenceperCapita-Targeted-Diesel.mat');

for ii=1:1330
[~,beta(ii,:),~,~,DA(ii),~,~,KP(ii,:),KV(ii),dV(ii,:),~,DAR(ii),w(ii)]=RetParameterPS(par(ii,:),XU,CF,4);
end