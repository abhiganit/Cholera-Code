clear;
clc;
X=struct('N',{'WaSH and Targeted Attacks ','WaSH, Food security, and Conflict ','WaSH, Food security, and Shelling ','WaSH, Food security, and Diesel ','Food security and Wheat ','WaSH and Rainfall ','WaSH, Food security, and Incidence per capita' });

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});

CF=[1;1];
RF=2;    

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
INDXS=[];
INDX=[3];
    load(['Fit-Vaccination-IncidenceperCapita' C(unique([INDXS])).N '-CF=1.mat'],'par')
parf=par;
        XU=zeros(1,28);
        XU([25:28])=1; 
        for ii=1:length(INDX)  
            XU([1:4]+4.*(INDX(ii)-1))=1;   
        end 
                  
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
          TestError(mm)=(OFuncProGA(part,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW))); 
        end

        SR=sortrows([partest TestError],length(partest(1,:))+1);
        temppar=SR(1,1:(end-1));

        [par,RSSv] =ProFittingGA(XU,CF,RF,temppar);  
        parf=[parf;par];
        save(['Fit-Vaccination-IncidenceperCapita' C(unique([INDX])).N '-CF=1.mat'],'par','RSSv','XU','CF','RF','X');
