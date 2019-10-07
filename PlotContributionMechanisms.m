% Read table of past fitsclose all;
close all;
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

    CInc=(beta(1).*squeeze(X(1,:,:))+beta(2).*squeeze(X(2,:,:))+beta(3).*squeeze(X(3,:,:))+beta(4).*squeeze(X(4,:,:))+beta(9).*squeeze(X(9,:,:))+beta(7).*squeeze(X(7,:,:))+beta(8).*squeeze(X(8,:,:)));
    CConflict=(beta(5).*squeeze(X(5,:,:))+beta(6).*squeeze(X(6,:,:))+beta(10).*squeeze(X(10,:,:))+beta(11).*squeeze(X(11,:,:))+beta(12).*squeeze(X(12,:,:))+beta(13).*squeeze(X(13,:,:)));
    Mt=sum(Yt);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

%Tot=repmat((sum(CInc,1)+sum(CRain,1)+sum(CConflict,1))./(Mt),3,1);
h=bar([1:149]+maxtau,([sum(CInc,1);sum(CConflict,1)])','stacked','LineStyle','none');
h(1).FaceColor = 'flat';
h(1).CData = [0.4 0.4 0.4];
h(2).FaceColor = 'flat';
h(2).CData = [0.9 0 0];
% The size to separate the weeks in the x-label

hold on;
%hmf=plot([1:149]+maxtau,NatIData(1+maxtau:end),'k-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
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
text(yh.Extent(1),55528,'A','Fontsize',32,'FontWeight','bold');
print(gcf,[pwd '\Figures\Figure1A.png'],'-dpng','-r600');
%% Waves
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.075105042016807,0.163120567375887,0.259453781512605,0.8217]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

Tot1=sum(sum(CInc(:,IW(1,1):IW(1,2)),1)+sum(CRain(:,IW(1,1):IW(1,2)),1)+sum(CConflict(:,IW(1,1):IW(1,2)),1));
Tot2=sum(sum(CInc(:,IW(2,1):IW(2,2)),1)+sum(CRain(:,IW(2,1):IW(2,2)),1)+sum(CConflict(:,IW(2,1):IW(2,2)),1));
Tot3=sum(sum(CInc(:,IW(3,1):IW(3,2)),1)+sum(CRain(:,IW(3,1):IW(3,2)),1)+sum(CConflict(:,IW(3,1):IW(3,2)),1));
Tot4=sum(sum(CInc(:,IW(4,1):IW(4,2)),1)+sum(CRain(:,IW(4,1):IW(4,2)),1)+sum(CConflict(:,IW(4,1):IW(4,2)),1));
h=barh([1:4],flip([sum(sum(CConflict(:,IW(1,1):IW(1,2)),1))./Tot1 sum(sum(CRain(:,IW(1,1):IW(1,2)),1))./Tot1;sum(sum(CConflict(:,IW(2,1):IW(2,2)),1))./Tot2 sum(sum(CRain(:,IW(2,1):IW(2,2)),1))./Tot2;sum(sum(CConflict(:,IW(3,1):IW(3,2)),1))./Tot3 sum(sum(CRain(:,IW(3,1):IW(3,2)),1))./Tot3;sum(sum(CConflict(:,IW(4,1):IW(4,2)),1))./Tot4 sum(sum(CRain(:,IW(4,1):IW(4,2)),1))./Tot4]),'stacked','LineStyle','none');

h(2).FaceColor = 'flat';
h(2).CData = [0 0.6 1];
h(1).FaceColor = 'flat';
h(1).CData = [0.9 0 0];
legend({'Conflict','Rainfall'},'Fontsize',18);
legend boxoff;
N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});

% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'YTick',[1:4],'XTick',[0:0.1:1],'Xminortick','on','Fontsize',16,'YTickLabel',{'Fourth wave','Third wave','Second wave','First wave'});

box off;
xtickangle(45);

xlabel({'Estimated contribtuion to suscpeted','cholera cases'},'Fontsize',18);
xlim([0 0.4]);
ylim([0.5 4.5]);
text(-0.290688259109312*max(xlim),4.443896424167694,'B','Fontsize',32,'FontWeight','bold');

%% Govnorate level various waves and the impact of rainfall
% 
Tot1=(sum(CInc(:,IW(1,1):IW(1,2)),2)+sum(CRain(:,IW(1,1):IW(1,2)),2)+sum(CConflict(:,IW(1,1):IW(1,2)),2));
Tot2=(sum(CInc(:,IW(2,1):IW(2,2)),2)+sum(CRain(:,IW(2,1):IW(2,2)),2)+sum(CConflict(:,IW(2,1):IW(2,2)),2));
Tot3=(sum(CInc(:,IW(3,1):IW(3,2)),2)+sum(CRain(:,IW(3,1):IW(3,2)),2)+sum(CConflict(:,IW(3,1):IW(3,2)),2));
Tot4=(sum(CInc(:,IW(4,1):IW(4,2)),2)+sum(CRain(:,IW(4,1):IW(4,2)),2)+sum(CConflict(:,IW(4,1):IW(4,2)),2));
% 
TT=[Tot1 Tot2 Tot3 Tot4];
R=[sum(CRain(:,IW(1,1):IW(1,2)),2)./Tot1 sum(CRain(:,IW(2,1):IW(2,2)),2)./Tot2 sum(CRain(:,IW(3,1):IW(3,2)),2)./Tot3 sum(CRain(:,IW(4,1):IW(4,2)),2)./Tot4];
CM=[sum(CConflict(:,IW(1,1):IW(1,2)),2)./Tot1 sum(CConflict(:,IW(2,1):IW(2,2)),2)./Tot2 sum(CConflict(:,IW(3,1):IW(3,2)),2)./Tot3 sum(CConflict(:,IW(4,1):IW(4,2)),2)./Tot4;];
for jj=1:4
    subplot('Position',[0.38,0.163120567375887+0.21*(jj-1),0.57,0.189]);
    yyaxis left
    h=bar([1:20],[CM(:,(5-jj)) R(:,(5-jj))],'stacked','LineStyle','none');
    h(2).FaceColor = 'flat';
    h(2).CData = [0 0.6 1];
    h(1).FaceColor = 'flat';
    h(1).CData = [0.9 0 0];
    box off;
    xlim([0.5 20.5]);
    ylim([0 1]);
    if(jj==1)        
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',[0:0.2:1],'Yminortick','on','Fontsize',16,'XTickLabel',{N(GNZI).G},'ycolor','k');
        xlabel('Govnorate','Fontsize',18)   
        xtickangle(45);
    else
        if(jj==2)            
            yh=ylabel({'Estimated contribtuion to suscpeted cholera cases'},'Fontsize',18);
            yh.Position=[-0.416870970915301,1.002674273628602,-1];
        end
        if(jj==4)
            text(-1,0.9591,'C','Fontsize',32','FontWeight','bold');
        end
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',[0:0.2:1],'Yminortick','on','Fontsize',16,'XTickLabel',{''},'ycolor','k');
    end
    yyaxis right
    semilogy([1:20],TT(:,5-jj),'-o','color',[0 0 0],'LineWidth',2,'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
    if(jj==1)        
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',10.^[0:6],'Fontsize',16,'XTickLabel',{N(GNZI).G},'ycolor',[0 0 0]);
        xlabel('Govnorate','Fontsize',18)   
        xtickangle(45);
    else        
        if(jj==2)  
            yh=ylabel({'Estimated suscpeted cholera cases'},'Fontsize',18);
            yh.Rotation=270;
            yh.Position=[22.13655849975374,3383877.445489365,-1];
        end
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',10.^[0:6],'Fontsize',16,'XTickLabel',{''},'ycolor',[0 0 0]);
    end
    box off;
    
    xlim([0.5 20.5]);
    ylim([1,10^6]);
end

print(gcf,[pwd '\Figures\Figure1B.png'],'-dpng','-r600');

