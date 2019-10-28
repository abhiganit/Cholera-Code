function [XUr,RSSr,CVEr,parr,kr,par] = ForwardSelection(XU,RSS,CVE,k,atest,PDS,pars)
%BACKWARDSSELECTION Takes the comblex model XU and determines if a simplier
%model is more suitble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify the model being fit    
        % XU(1)- beta_0
        % XU(2) - population density for the govenrate
        % XU(3)- Health facilities in the govnorates
        % XU(4) - Incidence last week
        % XU(5) - Inicedence in the other govnorates
        % XU(6) - Internally displaced people
        % XU(7) - Rebel control
        % XU(8) - Product of incidence and attacks 
        % XU(9) - Product of incidence and conflict 
        % XU(10) - Product of cumulative attacks incidence and rainfall 
        % XU(11)- cumulative attacks and Rainfall
        % XU(12) - Conflict, rianfall incidence
        % XU(13) - Attack, rianfall, incidence
% RSS - residual sum of squares of the model
% k - the number of paramters in the complex model
% atest - level of significance want to test
% PDS - percentage of data set to use in the fitting (0<=PDS<=1)
% pars - the starting point for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XUr - Specify the model being fit    
        % XU(1)- beta_0
        % XU(2) - population density for the govenrate
        % XU(3)- Health facilities in the govnorates
        % XU(4) - Incidence last week
        % XU(5) - Inicedence in the other govnorates
        % XU(6) - Internally displaced people
        % XU(7) - Rebel control
        % XU(8) - Product of incidence and attacks 
        % XU(9) - Product of incidence and conflict 
        % XU(10) - Product of cumulative attacks incidence and rainfall 
        % XU(11)- cumulative attacks and Rainfall
        % XU(12) - Conflict, rianfall incidence
        % XU(13) - Attack, rianfall, incidence
% RSSr - residual sum of squares of the model
% parr- the paramter of the reduced model
% kr - the number of paramters estimatedi in the reduced model

%% Start the backwards selection 
frindx=find(XU==0);
NMR=length(frindx);
XUm=repmat(XU,NMR,1);
par=zeros(NMR,length(pars(1,:)));
RSSv=zeros(NMR,1);
kr=zeros(NMR,1);
CVEv=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=1; % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVEv(ii)] = ProFittingGA(XUm(ii,:),PDS,pars);
    [kr(ii)]=RetParameterPS(par(ii,:),XUm(ii,:));
end
if(atest~=0)
     load('Yemen_Gov_Incidence.mat'); % Incidence data
    WI=IData'; % Transpose the data set such that the number of areas is the row
    maxtau=4;
    %Find areas where we have non-zero incidence over course of epidemic
    GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
     NWF=153;
%     NGS=floor(length(GNZI)*PDS);
%     Itemp=sum(WI(GNZI,:),2);
%     GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
%     % Find the top max
%     for ii=1:ceil(NGS/2)
%        f=find(Itemp==max(Itemp)); % Find the maximum
%        Itemp(f)=0; % set maximum to zero for it is no longer selected
%        GTF(ii)=f; % Record index
%     end
%     Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
%     % Find the minimum contributors
%     for ii=(ceil(NGS/2)+1):NGS
%        f=find(Itemp==min(Itemp)); % Select minimum
%        Itemp(f)=max(Itemp); % Set to maximum for not selected again
%        GTF(ii)=f; % Record index
%     end
%     GTF=sort(GTF)';
GTF=find(GNZI~=9);
    WI=WI(GNZI(GTF),(1+maxtau):NWF);
    N=length(WI(:));
    Fstatistic=((N-kr)./(kr-k)).*((RSS-RSSv)./(RSSv));
    CrC=1-fcdf(Fstatistic,kr-k,kr);
    f=find(CrC==min(CrC)); % choose minimium p-value as we want to increase the size of the model
    if(min(CrC)<=atest) % Larger model accpeted if 
        XUr=XUm(f,:);
        RSSr=RSSv(f);
        parr=par(f,:);      
        kr=kr(f);
        CVEr=CVEv(f);
    else
        XUr=[];
        RSSr=[];
        parr=[];
        kr=[];
        CVEr=[];
    end
else % Does the model comparison based on the cross validation error
    f=find(CVEv==min(CVEv)); % choose minimium p-value as we want to increase the size of the model
    if(min(CVEv)<CVE) % Larger model accpeted if 
        XUr=XUm(f,:);
        RSSr=RSSv(f);
        parr=par(f,:);      
        kr=kr(f);
        CVEr=CVEv(f);
    else
        XUr=[];
        RSSr=[];
        parr=[];
        kr=[];
        CVEr=[];
    end
end
    
end

