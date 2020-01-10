clear;
clc;
load('Fit-Vaccination-PercentData=80_2.mat');

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,0.8);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
ndata=WI(GNZI(GTF),(maxtau+1):NW);
ndata=length(ndata(:));
AIC=zeros(length(CF(1,:)),1);
fmin=find(CVE==min(CVE));
if(length(fmin)>1)
    for ii=1:length(AIC)
        [k]=RetParameterPS(par(ii,:),XU,CF(:,ii),RF);
        AIC(ii)= AICScore(k,ndata,RSSv(ii).*ndata);
    end
    fmin=find(AIC==min(AIC));
end
[k,beta]=RetParameterPS(par(fmin,:),XU,CF(:,fmin),RF);
XU=[double(beta(XU==1)>10^(-10)) 1 1];
X=struct('N',{'WaSH','WaSH and Targeted Attacks','WaSH and Conflict','WaSH and Shelling','WaSH and Diesel','Food security','Food security and Conflict','Food security and Shellings','Food security and Diesel','Food security and Wheat','WaSH and Rainfall (S.I.)','WaSH and Rainfall (I.P.C.)' });

[~,~,lbps,ubps,~,partemp] = BoundsFitting(XU,par(fmin,:),CF(:,fmin),[2;2]);


SS=10^4;
SN=SS+1;

CF=CF(:,fmin);
RF=[0 1 2 0 1 2 0 1 2;
    0 0 0 1 1 1 2 2 2];
RSSv=zeros(9,1);
CVE=zeros(9,1);
partest=zeros(SN,length([-32.*ones(1,length(XU)) -32 -32 -32 -6 -6 -6 -6 -32 -32 -32 log(0.9) -32 -32 -32 -32 -32]));
for tt=1:10
    lhs=lhsdesign(SS,length(lbps)); % Use latin-hypercube sampling
    partest2=[partemp;repmat(lbps,SS,1)+repmat(ubps-lbps,SS,1).*lhs];
    for mm=1:SN
        [partest(mm,:)] = ExpandPar(partest2(mm,:),XU,CF,[2;2],1);
    end
    TestError=zeros(SN,9);

    for jj=1:9
        parfor mm=1:SN
          [~,~,~,~,~,part] = BoundsFitting(XU,partest(mm,:),CF,RF(:,jj));
          TestError(mm,jj)=(OFuncProPS(part,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF(:,jj))); 
       end
    end


    parfor ii=1:9
       SR=sortrows([partest TestError(:,ii)],length(partest(1,:))+1);
       temppar=SR(1:50,1:(end-1));
        [par(ii,:),RSSv(ii),CVE(ii)] =ProFittingGA(XU,0.8,CF,RF(:,ii),temppar);
    end
    for kk=1:9
       for ii=1:9                   
           [~,~,~,~,~,part] = BoundsFitting(XU,par(ii,:),CF,RF(:,kk));
           Rtemp=(OFuncProPS(part,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF(:,kk)));
           if(Rtemp<RSSv(kk))
              RSSv(kk)=Rtemp; 
              par(kk,:)=par(ii,:);
              CVE(kk)=(OFuncProPS(part,CF,WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),XU,maxtau,RC(GNZI(GTCV)),WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),RF(:,kk)));
           end
       end
    end
    for ii=1:9
        [~,~,~,~,~,partemp(ii,:)] = BoundsFitting(XU,par(ii,:),CF,[2;2]);
    end
    save(['Fit-Vaccination-Rainfall-PercentData=80.mat'],'par','RSSv','CVE','XU','X','CF','RF');
end
 