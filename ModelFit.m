close all;

%% Load the data
%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,GNZI,maxtau] = LoadYemenData;
PDS=0.8;
atest=0.05;
%% Forward selection
load(['ForwardSelection-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
XU=XUv(end,:);
par=parv(end,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);

%% Run the projection

%% Run the logistic model with the data

[Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,GNZI));
Err=abs(Yt-WI(GNZI,1+maxtau:end));
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

NWF=floor(153*PDS);
plot([1+maxtau:NWF],sum(Yt(:,1:(NWF-maxtau))),'k',[1:NWF],sum(WI(GNZI,1:NWF)),'k*','LineWidth',2); hold on;
plot([NWF:153],sum(Yt(:,(NWF-maxtau):end)),'r',[NWF+1:153],sum(WI(GNZI,NWF+1:end)),'r*','LineWidth',2); hold on;

dW=4;
NW=153;
startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:5000:55000],'Yminortick','on','Fontsize',16,'XTickLabel',XTL);

    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
box off;
xtickangle(45);

xlabel('Date','Fontsize',18);
yh=ylabel({'Suspected cholera cases','(National)'},'Fontsize',18);
ylim([0 55005]);
xlim([1 153.5]);
X=struct('N',{'Bias','Population Density, ','Health Facilities, ','WASH and Incidence, ','Incidence and Attacks, ','Incidence and Conflict, ','WASH,Rainfall and Incidence, ','WASH and Rainfall, ','External Incidence, ','Attacks, ','Rebel Control, ','Conflict and Rainfall, ','Attack and Rainfall'});
title([X(XUv(end,:)==1).N ' (\alpha= ' num2str(atest) ')'],'Fontsize',20);

% create smaller axes in top right, and plot on it
axes('Position',[.45 .75 .5 .2])
box on
plot([1+maxtau:NWF],median(Err(:,1:(NWF-maxtau)),1),'k','LineWidth',2); hold on;
patch( [1+maxtau:NWF flip(1+maxtau:NWF)], [prctile(Err(:,1:(NWF-maxtau)),25) flip(prctile(Err(:,1:(NWF-maxtau)),75))],'k','Facealpha',0.3','LineStyle','none');

plot([NWF:153],median(Err(:,(NWF-maxtau):end),1),'r','LineWidth',2); hold on;
patch( [NWF:153 flip(NWF:153)], [prctile(Err(:,(NWF-maxtau):end),25) flip(prctile(Err(:,(NWF-maxtau):end),75))],'r','Facealpha',0.3','LineStyle','none');
% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'Yminortick','on','Fontsize',14,'XTickLabel',XTL);
box off;
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
box off;
xtickangle(45);

xlim([1 153]);
xlabel('Date','Fontsize',16);
ylabel({'Absolute error','(Governorate)'},'Fontsize',16);
