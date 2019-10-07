function [XUr,RSSr,parr,kr] = BackwardsSelection(XU,RSS,k,atest,PDS)
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
par=zeros(NMR,31);
RSSv=zeros(NMR,1);
kr=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=0; % Remove the covariate from the model
    [par(ii,:),~,RSSv(ii)] = ProFittingGA(XUm(ii,:),PDS,[],0,0,0,[-16.*ones(1,11) ones(1,5) 0 0 0 -16.*ones(1,8)]);
    [kr(ii)]=RetParameterPS(par(ii,:),XUm(ii,:));
end
load('Yemen_Gov_Incidence.mat'); % Incidence data
WI=IData'; % Transpose the data set such that the number of areas is the row

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
maxtau=4;
NWF=floor(153*PDS);
WI=WI(GNZI,(NWF+1):end);
N=length(WI(:));
Fstatistic=((N-k)./(k-kr)).*((RSSv-RSS)./(RSS));
CrC=1-fcdf(Fstatistic,kr-k,kr);
f=find(CrC==max(CrC)); % choose maximum p-value as we want to decrease the size of the model
if(max(CrC)>=atest) % smaller model accpeted if does not meet the criteria
    XUr=XUm(f,:);
    RSSr=RSSv(f);
    parr=par(f,:);  
else
    XUr=[];
    RSSr=[];
    parr=[];
    kr=[];
end
    
end

