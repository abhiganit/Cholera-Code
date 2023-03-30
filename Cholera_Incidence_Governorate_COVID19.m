close all;
[WI2,~,tA2,Rtv2,Temptv2,~,P2,RC,H2,WPINm2,FPINm2,~,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
% 
load('PopulationSize_Yemen.mat');
load('CholeraIncidence_COVID-19.mat');
CC19=cell2mat(CholeraCOVID19(:,2:end))';
CC19Date=datenum(CholeraCOVID19(:,1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Incidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


w=0.83/3;
h=0.85/4;
yp=flip(linspace(0.05,1,5))+0.01;
xp=linspace(0.05,1,4);
yp=yp(2:end);
xp=xp(1:end-1);
% 
XTL={datestr(datenum('March 30, 2020')+7.*[-18:6:55],'mmm dd, yyyy')};

[yp,xp]=meshgrid(yp,xp);
xp=xp(:);
yp=yp(:);
for ii=1:22
    dx1=linspace(0,max(CC19(ii,37:end)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_Incidence_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([-18:55],CC19(ii,37:end),'LineStyle','none','FaceColor',hex2rgb('#C60000')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=-18:55
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<CC19(ii,37+18+mm)))),dx1(dx1<CC19(ii,37+18+mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<CC19(ii,37+18+mm)))),dx2(dx2<CC19(ii,37+18+mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(55,max(1+CC19(ii,37:end))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[-18:6:55],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[-18:6:55],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([-18.5 55.5]);
        ylim([0 max(1+CC19(ii,37:end))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([-18:55],CC19(ii,37:end),'LineStyle','none','FaceColor',hex2rgb('#C60000')); hold on;   
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=-18:55
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<CC19(ii,37+18+mm)))),dx1(dx1<CC19(ii,37+18+mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<CC19(ii,37+18+mm)))),dx2(dx2<CC19(ii,37+18+mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(55,max(1+CC19(ii,37:end))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[-18:6:55],'XTickLabel',XTL,'Xminortick','on');
        xtickangle(90);
        xlabel('Week of report','Fontsize',20);
        xlim([-18.5 55.5]);
        ylim([0 max(1+CC19(ii,37:end))]);
        box off;
    end
    if(ii==10)       
        ylabel('Incidence per 10,000','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Incidence_2.png'],'-dpng','-r600');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Weekly conflict
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_COVID-19_Timeline.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,258);


XTL={datestr(datenum('April 6, 2020')+7.*[165:6:238]-7.*184,'mmm dd, yyyy')};

for ii=1:22
    dx1=linspace(0,max(Ctv(ii,165:238)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_conflict_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([165:238],Ctv(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#DE7A22')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Ctv(ii,mm)))),dx1(dx1<Ctv(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Ctv(ii,mm)))),dx2(dx2<Ctv(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(234,max(1+Ctv(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([164.5 238.5]);
        ylim([0 max(1+Ctv(ii,165:238))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([165:238],Ctv(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#DE7A22')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Ctv(ii,mm)))),dx1(dx1<Ctv(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Ctv(ii,mm)))),dx2(dx2<Ctv(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(234,max(1+Ctv(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        
        xlim([164.5 238.5]);
        ylim([0 max(1+Ctv(ii,165:238))]);
        box off;
    end
    if(ii==10)       
        ylabel('Number of conflict events','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Conflict_2.png'],'-dpng','-r600');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Shelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

load('Yemen_Air_Shelling_COVID-19.mat');
Mt=GLevelConflict(YASt,S,258);


XTL={datestr(datenum('April 6, 2020')+7.*[165:6:238]-7.*184,'mmm dd, yyyy')};

for ii=1:22
    dx1=linspace(0,max(Mt(ii,165:238)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_shelling_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([165:238],Mt(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#4C3F54')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Mt(ii,mm)))),dx1(dx1<Mt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Mt(ii,mm)))),dx2(dx2<Mt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(234,max(1+Mt(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([164.5 238.5]);
        ylim([0 max(1+Mt(ii,165:238))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([165:238],Mt(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#4C3F54')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Mt(ii,mm)))),dx1(dx1<Mt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Mt(ii,mm)))),dx2(dx2<Mt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(234,max(1+Mt(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        
        xlim([164.5 238.5]);
        ylim([0 max(1+Mt(ii,165:238))]);
        box off;
    end
    if(ii==7)       
        ylabel('Number of shelling/attacks','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Shelling_2.png'],'-dpng','-r600');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Diesle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Dieselt=DieselCOVID19;
load('Diesel_Gov_Yemen_COVID-19.mat');
Dieselt=Dieselt+min(Diesel(Diesel>0));

XTL={datestr(datenum('April 6, 2020')+7.*[165:6:238]-7.*184,'mmm dd, yyyy')};

for ii=1:22
    dx1=linspace(0,max(Dieselt(ii,165:238)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_Diesel_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([165:238],Dieselt(ii,165:238),'LineStyle','none','FaceColor',[153,52,4]./255); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Dieselt(ii,mm)))),dx1(dx1<Dieselt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Dieselt(ii,mm)))),dx2(dx2<Dieselt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(166,max(1+Dieselt(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','left');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([164.5 238.5]);
        ylim([0 max(1+Dieselt(ii,165:238))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([165:238],Dieselt(ii,165:238),'LineStyle','none','FaceColor',[153,52,4]./255); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Dieselt(ii,mm)))),dx1(dx1<Dieselt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Dieselt(ii,mm)))),dx2(dx2<Dieselt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(166,max(1+Dieselt(ii,165:238))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','left');
        
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        
        xlim([164.5 238.5]);
        ylim([0 max(1+Dieselt(ii,165:238))]);
        box off;
    end
    if(ii==10)       
        ylabel('Price of Diesel','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Diesel_2.png'],'-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Wheat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Wheatt=WheatCOVID19;
load('Wheat_Gov_Yemen_COVID-19.mat');
Wheatt=Wheatt+min(Wheat(Wheat>0));

XTL={datestr(datenum('April 6, 2020')+7.*[165:6:238]-7.*184,'mmm dd, yyyy')};

for ii=1:22
    dx1=linspace(0,max(Wheatt(ii,165:238)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_Wheat_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([165:238],Wheatt(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#FAAF08')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Wheatt(ii,mm)))),dx1(dx1<Wheatt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Wheatt(ii,mm)))),dx2(dx2<Wheatt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(166,max(1+Wheatt(ii,165:238))*1.05,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','left');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([164.5 238.5]);
        ylim([0 1.3.*max(1+Wheatt(ii,165:238))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([165:238],Wheatt(ii,165:238),'LineStyle','none','FaceColor',hex2rgb('#FAAF08')); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:238
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Wheatt(ii,mm)))),dx1(dx1<Wheatt(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Wheatt(ii,mm)))),dx2(dx2<Wheatt(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(166,max(1+Wheatt(ii,165:238))*1.05,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','left');
        
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:238],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        
        xlim([164.5 238.5]);
        ylim([0 1.3.*max(1+Wheatt(ii,165:238))]);
        box off;
    end
    if(ii==10)       
        ylabel('Price of Wheat','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Wheat_2.png'],'-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Rainfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
load('Rainfall_COVID-19.mat');
Rtv=RainCOVID19';


XTL={datestr(datenum('April 6, 2020')+7.*[165:6:221]-7.*184,'mmm dd, yyyy')};

for ii=1:22
    dx1=linspace(0,max(Rtv(ii,165:221)),41);
    dx2=dx1+(dx1(2)-dx1(1))./2;
    dx2=dx2(1:end-1);
    if(ii==1)||(ii==13)
        if(ii==13)
            print(gcf,['COVID-19_Rainfall_1.png'],'-dpng','-r600');
        end
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if(ii>=13)
        iin=ii-12;
    else
        iin=ii;
    end
    if(ii<21)
        subplot('Position',[xp(iin) yp(iin) w h])
        b=bar([165:221],Rtv(ii,165:221),'LineStyle','none','FaceColor',[5,112,176]./255); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:221
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Rtv(ii,mm)))),dx1(dx1<Rtv(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Rtv(ii,mm)))),dx2(dx2<Rtv(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(221,max(1+Rtv(ii,165:221))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        if(ii<19)
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:221],'XTickLabel',{},'Xminortick','on');
        else
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:221],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        end
        xlim([164.5 221.5]);
        ylim([0 max(1+Rtv(ii,165:221))]);
        box off;
    elseif(ii==22)        
        subplot('Position',[xp(9) yp(9) w h])
        b=bar([165:221],Rtv(ii,165:221),'LineStyle','none','FaceColor',[5,112,176]./255); hold on;      
        
        b.FaceColor = 'flat';
        for mm=1:19
            b.CData(mm,:) = [0.7 0.7 0.7];
        end
        if(RC(ii)==1)
            for mm=165:221
                 for jj=1:3
                    if((jj==1)||(jj==3))
                       scatter((0.3.*(jj-2)+ mm).*ones(size(dx1(dx1<Rtv(ii,mm)))),dx1(dx1<Rtv(ii,mm)),2,'k','filled');
                    else
                        scatter((0.3.*(jj-2)+ mm).*ones(size(dx2(dx2<Rtv(ii,mm)))),dx2(dx2<Rtv(ii,mm)),2,'k','filled');
                    end
                 end
            end
        end
        text(221,max(1+Rtv(ii,165:221))*0.95,S(ii).ADM1_EN,'Fontsize',18,'HorizontalAlignment','right');
        
            set(gca,'LineWidth',2,'tickdir','out','fontsize',20,'XTick',[165:6:221],'XTickLabel',XTL,'Xminortick','on');
            xtickangle(90);
            xlabel('Week of report','Fontsize',20);
        
        xlim([164.5 221.5]);
        ylim([0 max(1+Rtv(ii,165:221))]);
        box off;
    end
    if(ii==10)       
        ylabel('Rainfall (mm)','Fontsize',20); 
    end
end
print(gcf,['COVID-19_Rainfall_2.png'],'-dpng','-r600');



