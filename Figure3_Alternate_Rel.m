%% Inlcudes the effwects of conflict and shellngs on the diesel prices
% Read table of past fitsclose all;
close all;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);

[Yt,X]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
 CCR=cell(6,1);
 for mm=1:6
     tempmat=zeros(size(squeeze(X(1,:,:))));
    for ii=(maxtau*(mm-1)+1):(mm.*maxtau)
        tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000;
    end
    CCR{mm}=tempmat;
 end
%% Conflict indirect effect
load('DieselrepresentedthroughConflictShellings.mat','bd','XC','XS');
mmt=4;
tempmat=zeros(size(squeeze(X(1,:,:))));
tempmat2=zeros(size(squeeze(X(1,:,:))));
for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
    tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:)));
    tempmat2=tempmat2+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:)));
end
CCR{2}=CCR{2}+tempmat;
CCR{3}=CCR{3}+tempmat2;
CCR{4}=CCR{4}-tempmat-tempmat2;

IndW=[1 21; 22 74; 75 121; 122 149]; % Index of wave for the data used in the regression model
WW=zeros(4,length(GNZI),6);
for ww=1:4
    for mm=1:6
       WW(ww,:,mm) = (sum((CCR{mm}(:,IndW(ww,1):IndW(ww,2))),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2));
    end   
end

for mm=1:6
   CCR{mm}=CCR{mm}./MI; 
end

IndW=[1 21; 22 74; 75 121; 122 149]; % Index of wave for the data used in the regression model

%% Plot the data

ColorM=[[221,28,119]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        ]; 


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
Gintv=[ 1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18];
Gintv=[Gintv(:);19;20;21];
for nn=1:6
    figure('units','normalized','outerposition',[0 0 1 1]);

    for mm=1:3
        subplot('Position',[0.06,0.14+(3-mm)*0.29,0.935,0.27]);

        Gint=Gintv(mm+3.*(nn-1));
        b=bar([(1+maxtau):NW],[squeeze(CCR{1}(Gint,:)); squeeze(CCR{2}(Gint,:)); squeeze(CCR{3}(Gint,:)); squeeze(CCR{4}(Gint,:)); squeeze(CCR{5}(Gint,:)); squeeze(CCR{6}(Gint,:))]','Stacked','LineStyle','none');
        for ii=1:length(ColorM(:,1))
            b(ii).FaceColor = 'flat';
            b(ii).CData = ColorM(ii,:);
        end
        yh=ylabel('Relative contribution','Fontsize',18);
        xlim([0.5 NW+0.5]);
         ylim([0 8500]./7000*1);
        text(NW+0.5,1,XGL(Gint),'Fontsize',16,'HorizontalAlignment','right');
        set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',16,'XMinortick','off','YTick',[0:0.1:1]);
        box off;
        hold on;
        for ii=1:3
            plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,3500,1001),'k-.','LineWidth',2);        
        end
text(yh.Extent(1),max(ylim)*0.975,char(64+Gint),'Fontsize',32,'fontweight','bold');

        for wv=1:4
            for jj=1:length(squeeze(WW(1,1,:)))
                if(jj==1)
                    xstart=maxtau+IndW(wv,1);
                    xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*squeeze(WW(wv,Gint,jj));
                else                
                    xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj-1))));
                    xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:jj)));                
                end
                patch([xstart xstart xend xend],[max(ylim) 7250./7000*1 7250./7000*1 max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
                if(round(100.*squeeze(WW(wv,Gint,jj)))>=5)
                ht=text(mean([xstart xend]),mean([7250 8500]./7000*1),[num2str(round(100.*squeeze(WW(wv,Gint,jj)))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
                set(ht,'Rotation',90);
                end
            end
             xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj))));
             xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));   
             patch([xstart xstart xend xend],[max(ylim) 7250./7000*1 7250./7000*1 max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
             if(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))>=5)
             ht=text(mean([xstart xend]),mean([7250 8500]./7000*1),[num2str(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
             set(ht,'Rotation',90);
             end
         end
    end
    print(gcf,['Figure3SI_Relative' num2str(nn) '.png'],'-dpng','-r600');
end
nn=7;

    for mm=1:3
    figure('units','normalized','outerposition',[0 0 1 1]);
        subplot('Position',[0.06,0.14+(3-mm)*0.29,0.935,0.27]);

        Gint=Gintv(mm+3.*(nn-1));
        b=bar([(1+maxtau):NW],[squeeze(CCR{1}(Gint,:)); squeeze(CCR{2}(Gint,:)); squeeze(CCR{3}(Gint,:)); squeeze(CCR{4}(Gint,:)); squeeze(CCR{5}(Gint,:)); squeeze(CCR{6}(Gint,:))]','Stacked','LineStyle','none');
        for ii=1:length(ColorM(:,1))
            b(ii).FaceColor = 'flat';
            b(ii).CData = ColorM(ii,:);
        end
        yh=ylabel('Relative contribution','Fontsize',18);
        xlim([0.5 NW+0.5]);
         ylim([0 8500]./7000*1);
        text(NW+0.5,1,XGL(Gint),'Fontsize',16,'HorizontalAlignment','right');
        set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',16,'XMinortick','off','YTick',[0:0.1:1]);
        box off;
        hold on;
        for ii=1:3
            plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,3500,1001),'k-.','LineWidth',2);        
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
                patch([xstart xstart xend xend],[max(ylim) 7250./7000*1 7250./7000*1 max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
                if(round(100.*squeeze(WW(wv,Gint,jj)))>=5)
                ht=text(mean([xstart xend]),mean([7250 8500]./7000*1),[num2str(round(100.*squeeze(WW(wv,Gint,jj)))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
                set(ht,'Rotation',90);
                end
            end
             xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WW(wv,Gint,1:(jj))));
             xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));    
             patch([xstart xstart xend xend],[max(ylim) 7250./7000*1 7250./7000*1 max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
             if(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))>=5)
             ht=text(mean([xstart xend]),mean([7250 8500]./7000*1),[num2str(round(100.*(1-sum(squeeze(WW(wv,Gint,1:(jj))))))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
             set(ht,'Rotation',90);
             end
        end
         text(yh.Extent(1),max(ylim)*0.975,char(64+Gint),'Fontsize',32,'fontweight','bold');
set(gca,'linewidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'XMinortick','on','YTick',[0:0.1:1]);
xlabel('Week of report','Fontsize',18);
xtickangle(45);
box off;
print(gcf,['Figure3SI_Relative' num2str(nn) '_' num2str(mm) '.png'],'-dpng','-r600');
    end
