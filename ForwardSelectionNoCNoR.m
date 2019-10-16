function [XUr,RSSr,CVEr,parr,kr,par] = ForwardSelectionNoCNoR(XU,RSS,CVE,k,atest,PDS,Gov,pars)
%BACKWARDSSELECTION Takes the comblex model XU and determines if a simplier
%model is more suitble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU -  the X input that will be used in the fitting process. Will be zeros
 % or one (1 x 11)
    % XU(1)- beta_0
    % XU(2) - population density
    % XU(3) - number of health facilities 
    % XU(4) - Past incidence
    % XU(5) - Product of incidence and attacks
    % XU(6) - Product of incidence and conflict
    % XU(7) - Product of incidence and rainfall
    % XU(8) - Rainfall only        
    % XU(9) - Incidence in other govnorates
    % XU(10) - Attacks only
    % XU(11) - Rebel control
% RSS - residual sum of squares of the model
% k - the number of paramters in the complex model
% atest - level of significance want to test
% PDS - percentage of data set to use in the fitting (0<=PDS<=1)
% pars - the starting point for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XUr -  the X input that will be used in the fitting process. Will be zeros
 % or one (1 x 11)
    % XU(1)- beta_0
    % XU(2) - population density
    % XU(3) - number of health facilities 
    % XU(4) - Past incidence
    % XU(5) - Product of incidence and attacks
    % XU(6) - Product of incidence and conflict
    % XU(7) - Product of incidence and rainfall
    % XU(8) - Rainfall only        
    % XU(9) - Incidence in other govnorates
    % XU(10) - Attacks only
    % XU(11) - Rebel control
% RSSr - residual sum of squares of the model
% parr- the paramter of the reduced model
% kr - the number of paramters estimatedi in the reduced model

%% Start the backwards selection 

frindx=find(XU==0);
frindx=frindx(frindx<=5);
NMR=length(frindx);
XUm=repmat(XU,NMR,1);
par=zeros(NMR,length(pars));
RSSv=zeros(NMR,1);
kr=zeros(NMR,1);
CVEv=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=1; % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVEv(ii)] = ProFittingGA(XUm(ii,:),PDS,Gov,pars);    
    [kr(ii)]=RetParameterPS(par(ii,:),XUm(ii,:));
end
if(atest~=0)    
    N=floor(153*PDS);
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

