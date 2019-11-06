%% Produces the model prediction for national incidence and the cross validation governoertaes
close all;

%% Load the data
%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,IDPt,GNZI,maxtau] = LoadYemenData;
PDS=0.8;
atest=0.01;
%% Forward selection
load(['ModelSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
XU=XUv(end,:);
par=parv(end,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);

%% Determine the areas used in the fitting
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
GTF=sort(GTF)'; % Gov. to used in the fitting of the model. We sort to keep order consistent with GNZI

ErrCv=zeros(3,153-maxtau);

%% Run the model with all factors
[YtA,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,GNZI),IDPt(GNZI,:));
ErrFit=sum(abs(YtA(GTF,:)-WI(GNZI(GTF),1+maxtau:end)),1);
ErrCV(1,:)=sum(abs(YtA-WI(GNZI,1+maxtau:end)),1)-ErrFit;

%% Run Model with conflict and no rain
load(['ModelSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
XU=XUv(end-1,:);
par=parv(end-1,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
[YtC,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,GNZI),IDPt(GNZI,:));
ErrFit=sum(abs(YtC(GTF,:)-WI(GNZI(GTF),1+maxtau:end)),1);
ErrCV(2,:)=sum(abs(YtC-WI(GNZI,1+maxtau:end)),1)-ErrFit;

%% Run Model with rain and no conflict
load(['ModelSelectionNoRain-PercentDataSet=' num2str(PDS*100) '-alpha=' num2str(atest*100) '.mat']);
XU=XUv(end-2,:);
par=parv(end-2,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
[YtR,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,GNZI),IDPt(GNZI,:));
ErrFit=sum(abs(YtR(GTF,:)-WI(GNZI(GTF),1+maxtau:end)),1);
ErrCV(3,:)=sum(abs(YtR-WI(GNZI,1+maxtau:end)),1)-ErrFit;

ErrCV=ErrCV./(length(GNZI)-length(GTF));
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.16,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

NWF=153;
plot([1:NWF],sum(WI(GNZI,1:NWF)),'ko','LineWidth',2,'Markerfacecolor','k','Markersize',5); hold on;
plot([1+maxtau:NWF],sum(YtA(:,1:(NWF-maxtau))),'k','LineWidth',2);
plot([1+maxtau:NWF],sum(YtC(:,1:(NWF-maxtau))),'color',hex2rgb('#E81E25'),'LineWidth',2);
plot([1+maxtau:NWF],sum(YtR(:,1:(NWF-maxtau))),'color',hex2rgb('#6A92CC'),'LineWidth',2);

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

text(yh.Extent(1),0.99.*max(ylim),'A','Fontsize',32,'fontweight','bold');
% create smaller axes in top right, and plot on it
axes('Position',[.45 .75 .5 .2])
box on
plot([1+maxtau:NWF],ErrCV(1,:),'k','LineWidth',2); hold on;
plot([1+maxtau:NWF],ErrCV(2,:),'color',hex2rgb('#E81E25'),'LineWidth',2); hold on;
plot([1+maxtau:NWF],ErrCV(3,:),'color',hex2rgb('#6A92CC'),'LineWidth',2); hold on;

% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'Yminortick','on','Fontsize',14,'XTickLabel',XTL,'YTick',[0:100:1000]);
ylim([0 1000]);
box off;
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
box off;
xtickangle(45);

xlim([1 153]);
xlabel('Date','Fontsize',16);
ylabel({'Average Absolute error','(Governorate)'},'Fontsize',16);

%% Plot the cross validation for the governates
CVG=zeros(length(GNZI)-length(GTF),1);
cc=1;
for ii=1:length(GNZI)
    if(isempty(find(GNZI(ii)==GTF)))
        CVG(cc)=ii;
        cc=cc+1;
    end
end
figure('units','normalized','outerposition',[0 0 1 1]);
N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
for ii=1:2
    for jj=1:2
        subplot('Position',[0.0808+(0.94/2).*(jj-1),0.16+(2-ii).*(0.85/2),0.95.*0.897162184873949/2,0.95.*0.793313069908819/2]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
        plot([1:NWF],(WI(GNZI(CVG(jj+2.*(ii-1))),1:NWF)),'ko','LineWidth',2,'Markerfacecolor','k','Markersize',5); hold on;
        plot([1+maxtau:NWF],(YtA(CVG(jj+2.*(ii-1)),1:(NWF-maxtau))),'k','LineWidth',2);
        plot([1+maxtau:NWF],(YtC(CVG(jj+2.*(ii-1)),1:(NWF-maxtau))),'color',hex2rgb('#E81E25'),'LineWidth',2);
        plot([1+maxtau:NWF],(YtR(CVG(jj+2.*(ii-1)),1:(NWF-maxtau))),'color',hex2rgb('#6A92CC'),'LineWidth',2);
        title(N(GNZI(CVG(jj+2.*(ii-1)))).G,'Fontsize',18);
        dW=8;
        NW=153;
        startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
        XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
        % changing the aspects of the axis for the the current figure 
        if(ii==2)
            set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:500:5500],'Yminortick','on','Fontsize',16,'XTickLabel',XTL);
        xtickangle(45);

        xlabel('Date','Fontsize',18);
        else
            set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:500:5500],'Yminortick','on','Fontsize',16,'XTickLabel','');
        end
        if(jj==1)
        
        yh=ylabel({'Suspected cholera cases'},'Fontsize',18);
        end

            % Sets the y-axis to not have 10^n
            ax=gca; % finds the current axis
            ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        box off;
        ylim([0 5500]);
xlim([1 153.5]);
    end
end


