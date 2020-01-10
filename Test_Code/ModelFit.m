close all;

%% Load the data
%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,IDPt,GNZI,maxtau] = LoadYemenData;
PDS=0.8;
atest=0.05;
%% Forward selection
load(['ForwardSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
XU=XUv(end,:);
par=parv(end,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);

%% Determine the areas used in the fitting
NGS=floor(length(GNZI)*PDS);
Itemp=sum(WI(GNZI,:),2);
Itemp(9)=0; % set the one governorate of interest for analysis to zero so it is not included in the fitting but the cross validation
GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
% Find the top max
for ii=1:ceil(NGS/2)
   f=find(Itemp==max(Itemp)); % Find the maximum
   Itemp(f)=0; % set maximum to zero for it is no longer selected
   GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
Itemp(9)=max(Itemp); % set the one governorate of interest for analysis to maximum so it is not included in the fitting but the cross validation
% Find the minimum contributors
for ii=(ceil(NGS/2)+1):NGS
   f=find(Itemp==min(Itemp)); % Select minimum
   Itemp(f)=max(Itemp); % Set to maximum for not selected again
   GTF(ii)=f; % Record index
end
GTF=sort(GTF)'; % Gov. to used in the fitting of the model. We sort to keep order consistent with GNZI

%% Run the logistic model with the data

[Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),IDPt(GNZI,:));
ErrFit=sum(abs(Yt(GTF,:)-WI(GNZI(GTF),1+maxtau:end)),1);
ErrCV=sum(abs(Yt-WI(GNZI,1+maxtau:end)),1)-ErrFit;
ErrFit=ErrFit./length(GTF);
ErrCV=ErrCV./(length(GNZI)-length(GTF));
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

load('Yemen_Gov_Incidence.mat'); % Incidence data

WI=IData';

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
NWF=153;
plot([1+maxtau:NWF],sum(MI(:,1:(NWF-maxtau))),'k',[1:NWF],sum(WI(GNZI,1:NWF)),'k*','LineWidth',2); hold on;

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

% create smaller axes in top right, and plot on it
axes('Position',[.45 .75 .5 .2])
box on
plot([1+maxtau:NWF],ErrFit,'k','LineWidth',2); hold on;
plot([1+maxtau:NWF],ErrCV,'r','LineWidth',2); hold on;

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
ylabel({'Average Absolute error','(Governorate)'},'Fontsize',16);
