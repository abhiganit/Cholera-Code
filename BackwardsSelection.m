function [XUr,RSSr,CVEr,parr,kr,par] = BackwardsSelection(XU,RSS,CVE,k,atest,PDS,pars)
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
NMR=sum(XU);
frindx=find(XU==1);
XUm=repmat(XU,NMR,1);
par=zeros(NMR,length(pars(1,:)));
RSSv=zeros(NMR,1);
kr=zeros(NMR,1);
CVEv=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=0; % Remove the covariate from the model
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
    NGS=floor(length(GNZI)*PDS);
    Itemp=sum(WI(GNZI,:),2);
    GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
    % Find the top max
    for ii=1:ceil(NGS/2)
       f=find(Itemp==max(Itemp)); % Find the maximum
       Itemp(f)=0; % set maximum to zero for it is no longer selected
       GTF(ii)=f; % Record index
    end
    Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
    % Find the minimum contributors
    for ii=(ceil(NGS/2)+1):NGS
       f=find(Itemp==min(Itemp)); % Select minimum
       Itemp(f)=max(Itemp); % Set to maximum for not selected again
       GTF(ii)=f; % Record index
    end
    GTF=sort(GTF)';
    WI=WI(GNZI(GTF),(1+maxtau):NWF);
    N=length(WI(:));
    Fstatistic=((N-k)./(k-kr)).*((RSSv-RSS)./(RSS));
    CrC=1-fcdf(Fstatistic,kr-k,kr);
    f=find(CrC==max(CrC)); % choose maximum p-value as we want to decrease the size of the model
    if(max(CrC)>=atest) % smaller model accpeted if does not meet the criteria
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
else
    f=find(CVEv==max(CVEv)); % choose maximum p-value as we want to decrease the size of the model
    if(min(CVEv)<=CVE) % smaller model accpeted if does not meet the criteria
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

