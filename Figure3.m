%% Inlcudes the effwects of conflict and shellngs on the diesel prices
% Read table of past fitsclose all;
close all;
clear;

[IndW,WWRC,CCR,maxtau,GNZI,RC,~,~] = Contribution_to_Cholera_Incidence;
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
Gintv=[2 5 9 11 19 21];

CCRC=cell(7,1);
for ii=1:7
    Temp=[CCR{ii}];
    TempR=sum(Temp(RC(GNZI)==1,:),1);
    TempG=sum(Temp(RC(GNZI)==0,:),1);
    CCRC{ii}=[TempR;TempG];
end


figure('units','normalized','outerposition',[0 0 1 1]);
yyl=[1.2*10^4 1250];

yy0=[0.58 0.15];


shd=[1:8];
    labels = {'Target attacks','Weekly conflict','Shellings/attacks','Diesel','Wheat','Rainfall','Temprature'};
    
for mm=1:2
    subplot('Position',[0.07,yy0(mm),0.915,0.41]);

    Gint=Gintv(mm);
    b=bar([(1+maxtau):NW],[squeeze(CCRC{1}(mm,:)); squeeze(CCRC{2}(mm,:)); squeeze(CCRC{3}(mm,:)); squeeze(CCRC{4}(mm,:)); squeeze(CCRC{5}(mm,:)); squeeze(CCRC{6}(mm,:)) ; squeeze(CCRC{7}(mm,:))]','Stacked','LineStyle','none');
    for ii=1:length(ColorM(:,1))
        b(ii).FaceColor = 'flat';
        b(ii).CData = ColorM(ii,:);
        if(mm==1)
%             if(ii==2 || ii==3 || ii==4 ||ii==6)
                text(1,11100-1050.*shd(ii), labels{ii},'Fontsize',18,'Color',ColorM(ii,:));
%             end
        end
    end
    yh=ylabel('Suspected cases','Fontsize',18);
    xlim([0.5 NW+0.5]);
     ylim([0 8500]./7000.*(yyl(mm)));
     if(mm==1)
        text(NW+0.5,(yyl(mm)),'\it{\bf{Houthi}}','Fontsize',18,'HorizontalAlignment','right');
     else
         text(NW+0.5,(yyl(mm)),'\it{\bf{Government}}','Fontsize',18,'HorizontalAlignment','right');
     end
     if(mm==1)
        set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',18,'XMinortick','off','YminorTick','on','YTick',[0:2000:(yyl(mm))]);
     else
        set(gca,'linewidth',2,'tickdir','out','XTick','','XTickLabel','','Fontsize',18,'XMinortick','off','YminorTick','on','YTick',[0:250:(yyl(mm))]); 
     end
    box off;
    hold on;
    for ii=1:3
        plot((maxtau+mean([IndW(ii,2) IndW(ii+1,1)])).*ones(1001,1),linspace(0,(yyl(mm)),1001),'k-.','LineWidth',2);        
    end
    
    
    for wv=1:4
        for jj=1:length(squeeze(WWRC(1,1,:)))
            if(jj==1)
                xstart=maxtau+IndW(wv,1);
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*squeeze(WWRC(wv,mm,jj));
            else                
                xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,mm,1:(jj-1))));
                xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,mm,1:jj)));                
            end
            patch([xstart xstart xend xend],[max(ylim) 7250./7000*(yyl(mm)) 7250./7000*(yyl(mm)) max(ylim)] ,ColorM(jj,:),'Edgealpha',0)
            if(round(100.*squeeze(WWRC(wv,mm,jj)))>=5)
            ht=text(mean([xstart xend]),mean([7250 8500]./7000*(yyl(mm))),[num2str(round(100.*squeeze(WWRC(wv,mm,jj)))) '%'],'Fontsize',18,'HorizontalAlignment','center','color','w','Fontweight','bold');
            set(ht,'Rotation',90);
            end
        end
         xstart=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1)).*sum(squeeze(WWRC(wv,mm,1:(jj))));
         xend=maxtau+IndW(wv,1)+(IndW(wv,2)-IndW(wv,1));   
         patch([xstart xstart xend xend],[max(ylim) 7250./7000*(yyl(mm)) 7250./7000*(yyl(mm)) max(ylim)] ,[0.7 0.7 0.7],'Edgealpha',0)
         if(round(100.*(1-sum(squeeze(WWRC(wv,mm,1:(jj))))))>=5)
         ht=text(mean([xstart xend]),mean([7250 8500]./7000*(yyl(mm))),[num2str(round(100.*(1-sum(squeeze(WWRC(wv,mm,1:(jj))))))) '%'],'Fontsize',18,'HorizontalAlignment','center','color','w','Fontweight','bold');
         set(ht,'Rotation',90);
         end
    end
     text(yh.Extent(1)*0.9,max(ylim)*0.975,char(64+mm),'Fontsize',32,'fontweight','bold');
end

if(mm==1)
    set(gca,'linewidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',18,'XMinortick','on','YTick',[0:2000:(yyl(mm))]);
 else
    set(gca,'linewidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',18,'XMinortick','on','YTick',[0:250:(yyl(mm))]); 
 end
xlabel('Week of report','Fontsize',18);
xtickangle(45);
box off;
print(gcf,['Figure3'],'-depsc','-r600');
print(gcf,['Figure3'],'-dpng','-r600');