function [XUr,RSSr,parr,kr] = BackwardsSelection(XU,RSS,k,atest)
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
par=zeros(NMR,26);
RSSv=zeros(NMR,1);
for ii=1:NMR
    XUm(ii,frindx(ii))=0; % Remove the covariate from the model
    [par(ii,:),RSSv(ii)] = ProFittingGA(XUm(ii,:),[],0,0,0);
end
f=find(RSSv==min(RSSv));
RSSt=RSSv(f);
XUt=XUm(f,:);
[kr]=RetParameterGA(par(f,:),XUt);
load('Yemen_Gov_Incidence.mat'); % Incidence data
load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
WI=IData'; % Transpose the data set such that the number of areas is the row

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

NW=length(NatIData); % number of week want to go out to
Temp=WI(GNZI,:);
N=length(Temp(:))+NW-length(WI(1,:));
Fstatistic=((N-k)./(k-kr)).*((RSSt-RSS)./(RSS));
if(Fstatistic>finv(1-atest,k-kr,k))
    XUr=XUt;
    RSSr=RSSt;
    parr=par(f,:);    
else
    XUr=[];
    RSSr=[];
    parr=[];
    kr=[];
end
    
end

