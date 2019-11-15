% Read table of past fitsclose all;
close all;
load('Fit-Vaccination-PercentData=80.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,0.8);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
ndata=WI(GNZI(GTF),(maxtau+1):NW);
ndata=length(ndata(:));
AIC=zeros(length(CF(1,:)),1);
kv=zeros(length(CF(1,:)),1);
fmin=[1 2];%find(CVE==min(CVE));
if(length(fmin)>1)
    for ii=1:length(AIC)
        [kv(ii)]=RetParameterPS(par(ii,:),XU,CF(:,ii),RF);
        AIC(ii)= AICScore(kv(ii),ndata,RSSv(ii).*ndata);
    end
    fmin=find(AIC==min(AIC),1);
end
CF=CF(:,fmin);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,DBE,DAE,K,n,KP,KV,dV,r,KI,r0,rm]=RetParameterPS(par(fmin,:),XU,CF,RF);

[Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,tau,maxtau,CF,RC(GNZI),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),KP,V1(GNZI,1:length(WI(1,:))),V2(GNZI,1:length(WI(1,:))),KV,dV,r,KI,Rtv(GNZI,1:length(WI(1,:))),RF,r0,rm);
[YtNR,XNR]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,tau,maxtau,CF,RC(GNZI),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),KP,V1(GNZI,1:length(WI(1,:))),V2(GNZI,1:length(WI(1,:))),KV,dV,0,KI,Rtv(GNZI,1:length(WI(1,:))),RF,r0,rm);
dV1=ImpactAttack(V1(GNZI,1:length(WI(1,:)))-V2(GNZI,1:length(WI(1,:))),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,1:length(WI(1,:))),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
 

CIncWaSH=(1-EOVC).*beta(1).*(1+r.*repmat(RC(GNZI),1,length(PopS(GNZI,maxtau+1:end)))).*squeeze(X(1,:,:)).*PopS(GNZI,maxtau+1:end)./10000;
CCR=(1-EOVC).*beta(2).*(1+r.*repmat(RC(GNZI),1,length(PopS(GNZI,maxtau+1:end)))).*squeeze(X(2,:,:)).*PopS(GNZI,maxtau+1:end)./10000;
CIncFS=(1-EOVC).*beta(6).*(1+r.*repmat(RC(GNZI),1,length(PopS(GNZI,maxtau+1:end)))).*squeeze(X(6,:,:)).*PopS(GNZI,maxtau+1:end)./10000;
Rf=(1-EOVC).*beta(11).*(1+r.*repmat(RC(GNZI),1,length(PopS(GNZI,maxtau+1:end)))).*squeeze(X(11,:,:)).*PopS(GNZI,maxtau+1:end)./10000;
CRebC=(Yt-YtNR).*PopS(GNZI,maxtau+1:end)./10000;


IndW=[1 21; 22 74; 75 121; 122 149]; % Index of wave for the data used in the regression model
WW=zeros(4,length(GNZI),4);
for ww=1:4
    WW(ww,:,:)=[(sum(CIncFS(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2)) (sum(CIncWaSH(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2)) (sum(CCR(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2)) (sum(Rf(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2))];
end

%% Plot the data

ColorM=[hex2rgb('#2E4600'); % Food Security 
        0 0.6 0.6; % WaSH
        1 0 1;
        0 0 1]; % Targeted attacks


IW=7.*(([1; 22 ; 75 ; 122; 150]-1)+maxtau); % The 150 is the start of the week we do not have data for and we are subtracting a week for the index of the week as the index zero is Oct 3, 2016
IW=[IW(1) IW(2)-1 IW(2) IW(3)-1 IW(3) IW(4)-1 IW(4) IW(5)-1];
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+IW],'mmm.dd,yyyy');

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

XGL={S(GNZI).ADM1_EN};

dX=0.01;
Wid=0.2;
% Governorate


load('Yemen_Gov_Incidence.mat');
startDateofSim = datenum('10-03-2016');% Start date

IData=IData(:,GNZI)';

NW=length(IData(1,:));
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
Gintv=[2 5 9 13 19 21];
figure('units','normalized','outerposition',[0 0 1 1]);

for mm=1:3
    subplot('Position',[0.06,0.14+(3-mm)*0.29,0.935,0.27]);

    Gint=Gintv(mm);
    b=bar([(1+maxtau):NW],[(CIncFS(Gint,:)); CIncWaSH(Gint,:) ; CCR(Gint,:); Rf(Gint,:)]','Stacked','LineStyle','none');
    for ii=1:length(ColorM(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = ColorM(ii,:);
    end
    ylabel('Suspected cholera cases','Fontsize',18);
    xlim([0.5 NW+0.5]);
     ylim([0 8500]);
    text(NW+0.5,7000,XGL(Gint),'Fontsize',16,'HorizontalAlignment','right');
    set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',16,'XMinortick','off','YTick',[0:1000:7000]);
    box off;
    hold on;
    for ii=1:3
        plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,7000,1001),'k-.','LineWidth',2);        
    end
    
    
    for wv=1:4
        for jj=1:length(squeeze(WW(1,1,:)))
            if(jj==1)
                xstart=maxtau+IndW(wv,1);
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*squeeze(WW(wv,Gint,jj));
            else                
                xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj-1))));
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:jj)));                
            end
            patch([xstart xstart xend xend],[max(ylim) 7250 7250 max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
            if(round(100.*squeeze(WW(wv,Gint,jj)))>=5)
            ht=text(mean([xstart xend]),mean([7250 8500]),[num2str(round(100.*squeeze(WW(wv,Gint,jj)))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
            set(ht,'Rotation',90);
            end
        end
         xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj))));
         xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));   
         patch([xstart xstart xend xend],[max(ylim) 7250 7250 max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
         if(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))>=5)
         ht=text(mean([xstart xend]),mean([7250 8500]),[num2str(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
         set(ht,'Rotation',90);
         end
     end
end
figure('units','normalized','outerposition',[0 0 1 1]);
for mm=4:6
    subplot('Position',[0.06,0.14+(6-mm)*0.29,0.935,0.27]);

    Gint=Gintv(mm);
    b=bar([(1+maxtau):NW],[(CIncFS(Gint,:)); CIncWaSH(Gint,:) ; CCR(Gint,:); Rf(Gint,:)]','Stacked','LineStyle','none');
    for ii=1:length(ColorM(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = ColorM(ii,:);
    end
    ylabel('Suspected cholera cases','Fontsize',18);
    xlim([0.5 NW+0.5]);
    
     ylim([0 8500]);
    text(NW+0.5,7000,XGL(Gint),'Fontsize',16,'HorizontalAlignment','right');
    set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',16,'XMinortick','off','YTick',[0:1000:7000]);
    box off;
    hold on;
    for ii=1:3
        plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,7000,1001),'k-.','LineWidth',2);
    end
    
    
    for wv=1:4
        for jj=1:length(squeeze(WW(1,1,:)))
            if(jj==1)
                xstart=maxtau+IndW(wv,1);
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*squeeze(WW(wv,Gint,jj));
            else                
                xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj-1))));
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:jj)));                
            end
            patch([xstart xstart xend xend],[max(ylim) 7250 7250 max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
            if(round(100.*squeeze(WW(wv,Gint,jj)))>=5)
            ht=text(mean([xstart xend]),mean([7250 8500]),[num2str(round(100.*squeeze(WW(wv,Gint,jj)))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
            set(ht,'Rotation',90);
            end
        end
         xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj))));
         xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));   
         patch([xstart xstart xend xend],[max(ylim) 7250 7250 max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
         if(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))>=5)
         ht=text(mean([xstart xend]),mean([7250 8500]),[num2str(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
         set(ht,'Rotation',90);
         end
     end
end
set(gca,'linewidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'XMinortick','on','YTick',[0:1000:7000]);
xlabel('Week of report','Fontsize',18);
xtickangle(45);
box off;
% OLD VERSION that makes 4 different panels

% 
% % Governorate
% figure('units','normalized','outerposition',[0 0 1 1]);
% % First wave
% subplot('Position',[0.045,0.48,0.45,0.45]);
% 
% b=bar([1:length(GNZI)],squeeze(WW(1,:,:)),'LineStyle','none');
% set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
% xtickangle(45);
% box off;
% xlabel('Governorate','Fontsize',18);
% title([XTL(1,:) ' to ' XTL(2,:)],'Fontsize',18);
% for ii=1:length(ColorM(:,1))
%     b(ii).FaceColor = 'flat';
%     b(ii).CData = ColorM(ii,:);
% end
% ylim([0 1]);
% yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% ax=gca;
% tickl=ax.TickLength;
% text(yh.Extent(1),max(ylim)*1.05,'A','Fontsize',32,'FontWeight','bold');
% % Second wave
% subplot('Position',[0.545,0.48,0.45,0.45]);
% 
% b=bar([1:length(GNZI)],squeeze(WW(2,:,:)),'LineStyle','none');
% set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
% xtickangle(45);
% box off;
% xlabel('Governorate','Fontsize',18);
% title([XTL(3,:) ' to ' XTL(4,:)],'Fontsize',18);
% for ii=1:length(ColorM(:,1))
%     b(ii).FaceColor = 'flat';
%     b(ii).CData = ColorM(ii,:);
% end
% 
% ylim([0 1]);
% yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% 
% text(yh.Extent(1),max(ylim)*1.05,'B','Fontsize',32,'FontWeight','bold');
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% % Third wave
% subplot('Position',[0.045,0.48,0.45,0.45]);
% 
% b=bar([1:length(GNZI)],squeeze(WW(3,:,:)),'LineStyle','none');
% set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
% xtickangle(45);
% box off;
% xlabel('Governorate','Fontsize',18);
% title([XTL(5,:) ' to ' XTL(6,:)],'Fontsize',18);
% for ii=1:length(ColorM(:,1))
%     b(ii).FaceColor = 'flat';
%     b(ii).CData = ColorM(ii,:);
% end
% 
% ylim([0 1]);
% yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% 
% text(yh.Extent(1),max(ylim)*1.05,'C','Fontsize',32,'FontWeight','bold');
% % Fourth wave
% subplot('Position',[0.545,0.48,0.45,0.45]);
% 
% b=bar([1:length(GNZI)],squeeze(WW(4,:,:)),'LineStyle','none');
% set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
% xtickangle(45);
% box off;
% xlabel('Governorate','Fontsize',18);
% title([XTL(7,:) ' to ' XTL(8,:)],'Fontsize',18);
% for ii=1:length(ColorM(:,1))
%     b(ii).FaceColor = 'flat';
%     b(ii).CData = ColorM(ii,:);
% end
% 
% ylim([0 1]);
% yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% 
%     text(yh.Extent(1),max(ylim)*1.05,'D','Fontsize',32,'FontWeight','bold');
% % % National
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % subplot('Position',[0.045,0.48,0.955,0.45]);
% % 
% % 
% % b=bar([1:4],rand(4,3),'LineStyle','none');
% % ylim([0 1]);
% % set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',{'First','Second','Third','Fourth'},'Fontsize',16,'Yminortick','on','TickLength',tickl.*0.45/0.955);
% % 
% % box off;
% % xlabel('Epidemic wave','Fontsize',18);
% % ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% % for ii=1:length(ColorM(:,1))
% %     b(ii).FaceColor = 'flat';
% %     b(ii).CData = ColorM(ii,:);
% % end