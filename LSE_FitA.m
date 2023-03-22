pobj=parpool(32);
X=struct('N',{'WaSH and Targeted Attacks ','WaSH, Food security, and Conflict ','WaSH, Food security, and Shelling ','WaSH, Food security, and Diesel ','Food security and Wheat ','WaSH and Rainfall ','WaSH, Food Security and Temprature ','WaSH, Food security, and Incidence per capita' });

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});

[WI,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
% Loads the combinations of the different models
load('Combo.mat');

load(['LSE-Fit-Vaccination-IncidenceperCapita.mat'],'par');
par_base_t=par;

for nn=2:32
    % Construct the index used to determine the parameterizatrion
    indx=INC{nn};
    par_base=par_base_t;
    for kk=2:(nn-1)
        test_m=ismember(INC{kk},indx);
        if(sum(test_m)==length(INC{kk}))
            load(['LSE-Fit-Vaccination-IncidenceperCapita' C(INC{kk}).N '.mat'],'par');
            par_base=[par_base; par];
        end
    end
    
    XU=zeros(1,32);
    XU([21:32])=1; % Turn on rainfall, temprature and incidence
    for ii=1:length(indx)
        XU([1:4]+4.*(indx(ii)-1))=1;     
    end
    [lb,ub,~,~] = BoundsFitting(XU,[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 -32 log10(0.5) -32 -3],maxtau);
    % adjust bounds the searching bounds for the conflict function
    par_base_temp=zeros(size(par_base,1),length(lb));
    for kk=1:size(par_base,1)
        [~,~,~,par_base_temp(kk,:)] = BoundsFitting(XU,par_base(kk,:),maxtau);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Search for the parameterization to improve the optimziation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SR=[];
    for nss=1:100
        SS=5.*10^3;

        NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

        partest=zeros(SS,length([-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 -32 log10(0.5) -32 -3]));

        lhs=lhsdesign(SS,length(lb)); % Use latin-hypercube sampling
        if(nss==1)
            partest2=[par_base_temp; log10(repmat(10.^lb,SS,1)+repmat(10.^ub-10.^lb,SS,1).*lhs);repmat(lb,SS,1)+repmat(ub-lb,SS,1).*lhs];
            SS=2.*SS+size(par_base,1);
        else
            partest2=[log10(repmat(10.^lb,SS,1)+repmat(10.^ub-10.^lb,SS,1).*lhs);repmat(lb,SS,1)+repmat(ub-lb,SS,1).*lhs];
            SS=2.*SS;
        end
        
        parfor mm=1:SS
        [partest(mm,:)] = ExpandPar(partest2(mm,:),XU,maxtau);
        end
        TestError=zeros(SS,1);
        parfor mm=1:SS
          [~,~,~,part] = BoundsFitting(XU,partest(mm,:),maxtau);
          TestError(mm)=(OFuncProGA(part,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),Temptv(GNZI,1:NW),PopS(GNZI,1:NW),CI(GNZI,1:NW))); 
        end

        SRt=[sortrows([partest TestError],length(partest(1,:))+1)];
        SR=[SR; SRt(1:10,:)];
    end
    SR=sortrows(SR,length(partest(1,:))+1);
    temppar=SR(1:50,1:(end-1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Optimize based on the paramter stes temppar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    [par,RSSv] =ProFittingGA(XU,temppar);              
    save(['LSE-Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat'],'par','RSSv','XU','X');
end
delete pobj;