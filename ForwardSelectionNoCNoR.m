function [XUr,RSSr,CVEr,parr,kr,par] = ForwardSelectionNoCNoR(XU,RSS,k,atest,PDS,pars)
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
% No conflict
frindnr=find(frindx~=5);
frindx=frindx(frindnr);
frindnr=find(frindx~=6);
frindx=frindx(frindnr);
frindnr=find(frindx~=10);
frindx=frindx(frindnr);
frindnr=find(frindx~=11);
frindx=frindx(frindnr);
frindnr=find(frindx~=12);
frindx=frindx(frindnr);
frindnr=find(frindx~=13);
frindx=frindx(frindnr);
% No rainfall
frindnr=find(frindx~=7);
frindx=frindx(frindnr);
frindnr=find(frindx~=8);
frindx=frindx(frindnr);
NMR=length(frindx);
XUm=repmat(XU,NMR,1);
par=zeros(NMR,length(pars));
RSSv=zeros(NMR,1);
kr=zeros(NMR,1);
CVEv=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=1; % Remove the covariate from the model
    [par(ii,:),RSSv(ii),CVEv(ii)] = ProFittingGA(XUm(ii,:),PDS,pars);    
    [kr(ii)]=RetParameterPS(par(ii,:),XUm(ii,:));
end
load('Yemen_Gov_Incidence.mat'); % Incidence data
WI=IData'; % Transpose the data set such that the number of areas is the row
maxtau=4;
%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
NWF=floor(153*PDS);
WI=WI(GNZI,(1+maxtau):NWF);
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
    
end

