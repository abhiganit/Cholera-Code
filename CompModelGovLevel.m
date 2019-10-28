%% Compare Gov level incidence of the various models
close all;
clear;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,IDPt,GNZI,maxtau] = LoadYemenData;
PDS=0.8;
atest=0;
%% Forward selection
load(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
XU=XUr;
par=parr;

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
%% Run the projection

%% Run the logistic model with the data

[Yt,~]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),IDPt(GNZI,1:length(WI(1,:))));

% atest=0.05;
% %% Forward selection
% load(['ForwardSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
% XU=XUv(end,:);
% par=parv(end,:);
% 
% % Evaluate the number of paramters that are being used in the estimation 
% [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
% 
% %% Run the projection
% 

load(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
XU=XUv(end,:);
par=parv(end,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);

[YtNR,~]=LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),IDPt(GNZI,1:length(WI(1,:))));
NGS=floor(length(GNZI)*PDS);
Itemp=sum(WI(GNZI,:),2);
Itemp(RC(GNZI)==0)=0;
GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
% Find the top max
for ii=1:ceil(NGS/4)
   f=find(Itemp==max(Itemp)); % Find the maximum
   Itemp(f)=0; % set maximum to zero for it is no longer selected
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2);
Itemp(RC(GNZI)==1)=0;
% Find the top max
for ii=(ceil(NGS/4)+1):ceil(NGS/2)
   f=find(Itemp==max(Itemp)); % Find the maximum
   Itemp(f)=0; % set maximum to zero for it is no longer selected
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
Itemp(RC(GNZI)==0)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
% Find the minimum contributors
for ii=(ceil(NGS/2)+1):ceil(3.*NGS/4)
   f=find(Itemp==min(Itemp)); % Select minimum
   Itemp(f)=max(Itemp); % Set to maximum for not selected again
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
Itemp(RC(GNZI)==1)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
% Find the minimum contributors
for ii=(ceil(3.*NGS/4)+1):NGS
   f=find(Itemp==min(Itemp)); % Select minimum
   Itemp(f)=max(Itemp); % Set to maximum for not selected again
   GTF(ii)=f; % Record index
end
GTF=sort(GTF)'; % Gov. to used in the fitting of the model. We sort to keep order consistent with GNZI
load('Yemen_Gov_Incidence.mat'); % Incidence data

WI=IData';

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
MIr=(YtNR./(10000)).*PopS(GNZI,maxtau+1:end);
CompareGovFigureIncidence(MI,MIr,WI,GNZI,GTF)