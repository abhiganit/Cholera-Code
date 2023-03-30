% Read table of past fitsclose all;
close all;
Cov2Inc=[2 3 4 6];

[IndW,~,CCR,maxtau,GNZI,RC,WWRC] = Contribution_to_Cholera_Incidence;
%% Plot the data
ColorM=[[152,78,163]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        [247,129,191]./255;];   % Temprature


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
Gintv=[ 1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21];
Gintv=Gintv(:);

CVName={'Target attacks','Weekly conflict','Shellings/attacks','Diesel','Wheat','Rainfall','Temprature'};

shd=[0 0 1 2 2 3 4];
for Gint=1:length(GNZI)
    figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

    subplot('Position',[0.07,0.225,0.92,0.755]);

    b=bar([(1+maxtau):NW],[squeeze(CCR{1}(Gint,:)); squeeze(CCR{2}(Gint,:)); squeeze(CCR{3}(Gint,:)); squeeze(CCR{4}(Gint,:)); squeeze(CCR{5}(Gint,:)); squeeze(CCR{6}(Gint,:)); squeeze(CCR{7}(Gint,:))]','Stacked','LineStyle','none');
    for ii=1:length(ColorM(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = ColorM(ii,:);
    end
    myX=max(ylim).*1.1;
    yh=ylabel('Suspected cholera cases','Fontsize',18);
    xlim([0.5 NW+0.5]);
    ylim([0 8500]./7000*myX);
    text(1.15,myX,['\bf{\it{' XGL(Gint) '}}'],'Fontsize',16,'HorizontalAlignment','left');
    for ii=1:length(CVName)
        text(1.15,myX.*(1-0.06.*ii),CVName{ii},'Fontsize',16,'HorizontalAlignment','left','color',ColorM(ii,:));
    end
        
    set(gca,'linewidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'XMinortick','on');%,'YTick',[0:500:3500]);
    xlabel('Week of report','Fontsize',18);
    xtickangle(45);
    box off;
    hold on;
    for ii=1:3
        plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,myX,1001),'k-.','LineWidth',2);        
    end
%     text(yh.Extent(1),max(ylim)*0.975,char(64+Gint),'Fontsize',32,'fontweight','bold');

    for wv=1:4
        for jj=1:length(squeeze(WWRC(1,1,:)))
            if(jj==1)
                xstart=maxtau+IndW(wv,1);
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*squeeze(WWRC(wv,Gint,jj));
            else                
                xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,Gint,1:(jj-1))));
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,Gint,1:jj)));                
            end
            patch([xstart xstart xend xend],[max(ylim) 7250./7000*myX 7250./7000*myX max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
            if(round(100.*squeeze(WWRC(wv,Gint,jj)))>=5)
            ht=text(mean([xstart xend]),mean([7250 8500]./7000*myX),[num2str(round(100.*squeeze(WWRC(wv,Gint,jj)))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
            set(ht,'Rotation',90);
            end
        end
         xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,Gint,1:(jj))));
         xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));   
         patch([xstart xstart xend xend],[max(ylim) 7250./7000*myX 7250./7000*myX max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
         if(round(100.*(1-sum(squeeze(WWRC(wv,Gint,1:(jj))))))>=5)
         ht=text(mean([xstart xend]),mean([7250 8500]./7000*myX),[num2str(round(100.*(1-sum(squeeze(WWRC(wv,Gint,1:(jj))))))) '%'],'Fontsize',12,'HorizontalAlignment','center','color','w','Fontweight','bold');
         set(ht,'Rotation',90);
         end
    end
    
    print(gcf,['Gov_Contribution_' XGL{Gint} '.png'],'-dpng','-r600');
end
