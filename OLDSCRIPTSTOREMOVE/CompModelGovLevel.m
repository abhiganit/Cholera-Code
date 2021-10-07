%% Compare Gov level incidence of the various models
close all;
clear;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
PDS=0.8;
atest=0;
%% Forward selection
load('TestFit-Vaccination-PercentDataVARY.mat')
CF=[1;1];
RIF=zeros(4,1);
RF=0;
 [~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,mln,a,KV,dV]=RetParameterPS(par(1,:),XU,CF);

%% Run the projection

%% Run the logistic model with the data

[Yt,~]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a,V1(GNZI,1:length(WI(1,:))),V2(GNZI,1:length(WI(1,:))),KV,dV);
dI=WI(GNZI,(1+maxtau):end)-Yt;
temp=0;
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

beta(1)=0;
[YtNR,~]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a,V1(GNZI,1:length(WI(1,:))),V2(GNZI,1:length(WI(1,:))),KV,dV);

[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,PDS);
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