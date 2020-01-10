
pobj=parpool(20);
X=struct('N',{'WaSH and Targeted Attacks ','WaSH, Food security, and Conflict ','WaSH, Food security, and Shelling ','WaSH, Food security, and Diesel ','Food security and Wheat ','WaSH and Rainfall ','WaSH, Food security, and Incidence per capita' });

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});

    CF=[2 ;2];
    RF=1;    
    
    [WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
    load('Combo.mat');
    for nn=1:16
            indx=INC{nn};
            XU=zeros(1,28);
            XU([25:28])=1;
            for ii=1:length(indx)
                XU([1:4]+4.*(indx(ii)-1))=1;     
            end
            [~,~,lbps,ubps,~,~] = BoundsFitting(XU,[-32.*ones(1,length(XU)) log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32],CF,maxtau);
            % adjust bounds the searching bounds for the conflict function

            SS=5.*10^4;
            NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

            partest=zeros(SS,length([-32.*ones(1,length(XU)) log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32]));

            lhs=lhsdesign(SS,length(lbps)); % Use latin-hypercube sampling
            partest2=[repmat(lbps,SS,1)+repmat(ubps-lbps,SS,1).*lhs];
            parfor mm=1:SS
            [partest(mm,:)] = ExpandPar(partest2(mm,:),XU,CF,maxtau,1);
            end
            TestError=zeros(SS,1);
            parfor mm=1:SS
              [~,~,~,~,~,part] = BoundsFitting(XU,partest(mm,:),CF,maxtau);
              TestError(mm)=(OFuncProPS(part,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW))); 
            end

            SR=sortrows([partest TestError],length(partest(1,:))+1);
            temppar=SR(1:50,1:(end-1));

            [par,RSSv] =ProFittingGA(XU,CF,RF,temppar);              
            save(['Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat'],'par','RSSv','XU','CF','RF','X');
            
    end
            
            delete pobj;
 