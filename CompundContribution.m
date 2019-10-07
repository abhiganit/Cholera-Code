%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,GNZI,maxtau] = LoadYemenData;

    %% Forward selection
    load('ForwardSelection-PercentDataSet=80-alpha=1.mat');
    XU=XUv(end,:);
    par=parv(end,:);
    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
    %% Run the projection
    
    %% Run the logistic model with the data

    [Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI),WPIN(GNZI),Mt(GNZI,GNZI));

    CConflict=(beta(5).*squeeze(X(5,:,:))+beta(6).*squeeze(X(6,:,:))+beta(10).*squeeze(X(10,:,:))+beta(11).*squeeze(X(11,:,:)));
    CInc=(beta(2).*squeeze(X(2,:,:))+beta(3).*squeeze(X(3,:,:))+beta(4).*squeeze(X(4,:,:))+beta(9).*squeeze(X(9,:,:)));
    CRain=(beta(7).*squeeze(X(7,:,:))+beta(8).*squeeze(X(8,:,:)));
    CTotal=Yt;
    pr=zeros(20,153);
    pg=zeros(20,153);
    for ii=(1+maxtau):153
       pr(:,ii)=(CRain(:,ii-maxtau)+(beta(2).*pr(:,ii-tau(1)).*squeeze(X(2,:,ii-maxtau))'+beta(3).*pr(:,ii-tau(2)).*squeeze(X(3,:,ii-maxtau))'+beta(4).*pr(:,ii-tau(3)).*squeeze(X(4,:,ii-maxtau))'+beta(9).*pr(:,ii-tau(8)).*squeeze(X(9,:,ii-maxtau))'))./Yt(:,ii-maxtau);
       pg(:,ii)=(CConflict(:,ii-maxtau)+(beta(2).*pg(:,ii-tau(1)).*squeeze(X(2,:,ii-maxtau))'+beta(3).*pg(:,ii-tau(2)).*squeeze(X(3,:,ii-maxtau))'+beta(4).*pg(:,ii-tau(3)).*squeeze(X(4,:,ii-maxtau))'+beta(9).*pg(:,ii-tau(8)).*squeeze(X(9,:,ii-maxtau))'))./Yt(:,ii-maxtau);
       ff=find(Yt(:,ii-maxtau)==0);
       pr(ff,ii)=0;
       pg(ff,:)=0;
    end
    ContC=sum(Yt.*pg(:,1+maxtau:end));
    %ContR=sum(Yt.*pr(:,1+maxtau:end));
    ContO=sum(Yt)-ContC;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

%Tot=repmat((sum(CInc,1)+sum(CRain,1)+sum(CConflict,1))./(Mt),3,1);
h=bar([1:149]+maxtau,([ContO; ContC])','stacked','LineStyle','none');
h(1).FaceColor = 'flat';
h(1).CData = [0.4 0.4 0.4];
h(2).FaceColor = 'flat';
h(2).CData = [0.9 0 0];
% The size to separate the weeks in the x-label

hold on;
dW=4;
NW=length(WI(1,:));
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
yh=ylabel('Suspected cholera cases','Fontsize',18);
ylim([0 55005]);
xlim([1 153.5]);
IW=[1 21; 22 74; 75 116; 117 149];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',18);

legend([h],{'Other factors','Conflict'},'Fontsize',18,'location','northwest');

legend boxoff;

    test=0;