%% Runs the Fitting for the models
clear;
clc;
X=struct('N',{'WaSH and Targeted Attacks ','WaSH, Food security, and Conflict ','WaSH, Food security, and Shelling ','WaSH, Food security, and Diesel ','Food security and Wheat ','WaSH and Rainfall ','WaSH, Food security, and Incidence per capita' });

C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});

% Functions to be used for the conflict
CF=[2;2];
% Function to be sued for rainfall
RF=2;    

% Load the data to be used in the model and the fitting
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

% Index of the previous model to use as starting point
INDXS=[2 3 5 6];
% model to be fit
INDX=[1 2 3 5 6];
    load(['Fit-Vaccination-IncidenceperCapita' C(unique([INDXS])).N '.mat'],'par')
parf=par;

% set up the structure of the model with the indicators
        XU=zeros(1,28);
        XU([25:28])=1; 
        for ii=1:length(INDX)  
            XU([1:4]+4.*(INDX(ii)-1))=1;   
        end 
% The bounds to be used during the fitting process                  
        [~,~,lbps,ubps,~,~] = BoundsFitting(XU,[-32.*ones(1,length(XU))  -32 0 -1 0 -1 -32 -32 -32 -32 log(0.9) -32],CF,maxtau);
        % adjust bounds the searching bounds for the conflict function

        % Construct the starting points for the optimization algorithm
        % Number of points that will be used
        SS=5.*10^4;
        if(~isempty(parf))
            SN=SS+length(parf(:,1));
        else
            SN=SS;
        end
        
        % Number of weeks that will be used in the fitting
        NW=153; 

        % Allocate memory for the paramters
        partest=zeros(SS,length([-32.*ones(1,length(XU))  -32 0 -1 0 -1 -32 -32 -32 -32 log(0.9) -32]));
         
        %Sample space
        lhs=lhsdesign(SS,length(lbps)); % Use latin-hypercube sampling
        partest2=[repmat(lbps,SS,1)+repmat(ubps-lbps,SS,1).*lhs];
        % Construct the parameters expanding for the sampled bounds
        parfor mm=1:SS
            [partest(mm,:)] = ExpandPar(partest2(mm,:),XU,CF,maxtau,1);
        end
        % the points to be used in the fitting
        partest=[parf;partest];
        TestError=zeros(SN,1);
        % Allocate memory for the objective function
        parfor mm=1:SN
            % Returns parameters asscoaited with the specified model
          [~,~,~,~,~,part] = BoundsFitting(XU,partest(mm,:),CF,maxtau);
          % Computes the value of the objective function
          TestError(mm)=(OFuncProGA(part,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW))); 
        end

        SR=sortrows([partest TestError],length(partest(1,:))+1);
        temppar=SR(1,1:(end-1));

        [par,RSSv] =ProFittingGA(XU,CF,RF,temppar);  
        parf=[parf;par];
        save(['Fit-Vaccination-IncidenceperCapita' C(unique([INDX])).N '.mat'],'par','RSSv','XU','CF','RF','X');
