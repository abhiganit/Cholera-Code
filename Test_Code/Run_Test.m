clear;
clc;
X=struct('N',{'WaSH and Targeted Attacks ','WaSH, Food security, and Conflict ','WaSH, Food security, and Shelling ','WaSH, Food security, and Diesel ','Food security and Wheat ','WaSH and Rainfall ','WaSH, Food security, and Incidence per capita' });

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});

CF=[1 ;1];
RF=1;    

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,0.8);

% Cross validation Error Selected Model pair

load(['Fit-Vaccination-PercentData=80-IncidenceperCapita-Diesel-Rain.mat'],'par','RSSv','CVE','XU','CF','RF','X');

indx=[1 2 3 5];
parf=par;

for yy=1:4
    for ss=yy:4
        XU=zeros(1,28);
        XU([25:28])=1;
        XU([1:4]+4.*(4-1))=1;            
        XU([1:4]+4.*(6-1))=1;
                  
        XU([1:4]+4.*(indx(yy)-1))=1;                  
        XU([1:4]+4.*(indx(ss)-1))=1;
        [~,~,lbps,ubps,~,~] = BoundsFitting(XU,[-32.*ones(1,length(XU))  -32 0 -1 0 -1 -32 -32 -32 -32 log(0.9) -32],CF,maxtau);
        % adjust bounds the searching bounds for the conflict function

        SS=5.*10^4;
        if(~isempty(parf))
            SN=SS+length(parf(:,1));
        else
            SN=SS;
        end
        NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

        partest=zeros(SS,length([-32.*ones(1,length(XU))  -32 0 -1 0 -1 -32 -32 -32 -32 log(0.9) -32]));

        lhs=lhsdesign(SS,length(lbps)); % Use latin-hypercube sampling
        partest2=[repmat(lbps,SS,1)+repmat(ubps-lbps,SS,1).*lhs];
        parfor mm=1:SS
        [partest(mm,:)] = ExpandPar(partest2(mm,:),XU,CF,maxtau,1);
        end
        partest=[parf;partest];
        TestError=zeros(SN,1);
        parfor mm=1:SN
          [~,~,~,~,~,part] = BoundsFitting(XU,partest(mm,:),CF,maxtau);
          TestError(mm)=(OFuncProPS(part,CF,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF,PopS(GNZI(GTF),1:NW),CI(GNZI(GTF),1:NW))); 
        end

        SR=sortrows([partest TestError],length(partest(1,:))+1);
        temppar=SR(1,1:(end-1));

        [par,RSSv,CVE] =ProFittingGA(XU,0.8,CF,RF,temppar);  
        save(['Fit-Vaccination-PercentData=80-IncidenceperCapita-Diesel-Rain' C(unique([indx(ss) indx(yy)])).N '.mat'],'par','RSSv','CVE','XU','CF','RF','X');
    end
end