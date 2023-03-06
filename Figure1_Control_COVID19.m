close all;
clear;
[WI2,~,tA2,Rtv2,Temptv2,~,P2,RC2,H2,WPINm2,FPINm2,~,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
WI2=WI2(GNZI,:);
tA2=tA2(GNZI,:);
Rtv2=Rtv2(GNZI,:);
P2=P2(GNZI,:);
RC=RC2;
RC2=RC2(GNZI);
H2=H2(GNZI);
WPINm2=WPINm2(GNZI,:);
FPINm2=FPINm2(GNZI,:);
Wheatt2=Wheatt2(GNZI,:);


% Need to increase the diesel price as we transformed it by subtracting the
% minimum
Dieselt2=DieselCOVID19;
Dieselt2=Dieselt2(GNZI,:);

load('Diesel_Gov_Yemen_COVID-19.mat');
Dieselt2=Dieselt2+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt2=Wheatt2+min(Wheat(Wheat>0));

% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


load('Conflict_COVID-19_Timeline.mat'); % Load the conflict in the area for the projection
Ctv2=GLevelConflict(ProC,S,258);
load('Yemen_Air_Shelling_COVID-19.mat');
Mt2=GLevelConflict(YASt,S,258); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

Ctv2=Ctv2(GNZI,:);
Mt2=Mt2(GNZI,:);

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


figure('units','normalized','outerposition',[0 0 1 1]);
load('PopulationSize_Yemen.mat');
load('CholeraIncidence_COVID-19.mat');
CC19=cell2mat(CholeraCOVID19(:,2:end))';
CC19Date=datenum(CholeraCOVID19(:,1));
CC19=CC19(GNZI,:);

AP=AP(GNZI,:);
PopH=sum(AP(RC2==1,4));
PopG=sum(AP(RC2==0,4));
TestH=10000.*sum(CC19(RC2==1,:),1)./PopH;
TestG=10000.*sum(CC19(RC2==0,:),1)./PopG;


subplot('Position',[0.047268907563025,0.624113475177305,0.942752100840336,0.365886524822696]); 
plot([1:length(TestH)],TestH,'color',hex2rgb('#C60000'),'LineWidth',2); hold on;
plot([1:length(TestH)],TestG,'color',hex2rgb('#C60000'),'LineWidth',2,'LineStyle','-.');


plot(56.*ones(101,1),linspace(0,15,101),'k:','LineWidth',2);
startDateofSim = datenum('03-18-2019');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(length(TestH)-1)]],'mm/dd/yy');
xlim([0.5 length(TestH)+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:length(TestH)],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',18);
box off;
yh=ylabel({'Incidence per 10,000'},'Fontsize',18);
legend('Houthi control','Government control')
legend boxoff;

subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=10000.*mean(CC19(:,1:55)./repmat(AP(:,4),1,55),2);
MARo=10000.*mean(CC19(:,56:end)./repmat(AP(:,4),1,length(CC19(1,56:end))),2);
MART=MAR;
MAR=MAR./max(max(MAR),max(MARo));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   

box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(max(MART),max(MARo)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#C60000'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Incidence per 10,000','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

text(41.45,mean(ylim),'March 18, 2019 to April 5, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)


subplot('Position',[0.565,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=10000.*mean(CC19(:,56:end)./repmat(AP(:,4),1,length(CC19(1,56:end))),2);
MARo=10000.*mean(CC19(:,1:55)./repmat(AP(:,4),1,55),2);
MART=MAR;
MAR=MAR./max(max(MAR),max(MARo));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(max(MART),max(MARo)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#C60000'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Incidence per 10,000','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;


text(41.45,mean(ylim),'April 6,2020 to April 25, 2021','HorizontalAlignment','center','Rotation',90,'FontSize',20)

print(gcf,['Figure1_Incidence_COVID-19.png'],'-dpng','-r600');
%% Rainfall
load('Rainfall_COVID-19.mat');
Rtv2=RainCOVID19';
Rtv2=Rtv2(GNZI,:);
TestH=mean(Rtv2(RC2==1,:),1);
TestG=mean(Rtv2(RC2==0,:),1);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.047268907563025,0.624113475177305,0.942752100840336,0.365886524822696]); 

plot([146:221],TestH(146:221),'color',[5 112 176]./255,'LineWidth',2); hold on; 
plot([146:221],TestG(146:221),'color',[5 112 176]./255,'LineWidth',2,'LineStyle','-.'); hold on; 

plot(184.*ones(101,1),linspace(0,62,101),'k:','LineWidth',2);

startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(258-1)]],'mm/dd/yy');
xlim([146-0.5 221+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:258],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',18);
box off;
ylabel('Rainfall (mm)','Fontsize',18);
ylim([0 62])

legend('Houthi control','Government control')
legend boxoff;

subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Rtv2(:,[146:183]),2);
MARo=mean(Rtv2(:,[184:221]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[5 112 176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[5 112 176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii),2)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[5 112 176]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),2)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Rainfall (mm)','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

text(41.45,mean(ylim),'July 15, 2019 to April 5, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)


subplot('Position',[0.565,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Rtv2(:,[184:221]),2);
MARo=mean(Rtv2(:,[146:183]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[5 112 176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[5 112 176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii),2)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[5 112 176]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),2)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Rainfall (mm)','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;


text(41.45,mean(ylim),'April 6,2020 to Dec. 27, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)

print(gcf,['Figure1_Rain-COVID-19.png'],'-dpng','-r600');

%% Conflict


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Ctv2(RC2==1,:),1);
TestG=mean(Ctv2(RC2==0,:),1);


figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.047268907563025,0.624113475177305,0.942752100840336,0.365886524822696]); 

plot([109:258],TestH(109:258),'color',hex2rgb('#DE7A22'),'LineWidth',2); hold on; 
plot([109:258],TestG(109:258),'color',hex2rgb('#DE7A22'),'LineWidth',2,'LineStyle','-.'); hold on; 

plot(184.*ones(101,1),linspace(0,25,101),'k:','LineWidth',2);

startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(258-1)]],'mm/dd/yy');
xlim([109-0.5 258+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:258],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',18);
box off;
ylabel('Number of conflict events','Fontsize',18);

legend('Houthi control','Government control')
legend boxoff;
text(-13.9,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');


subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Ctv2(:,[109:183]),2);
MARo=sum(Ctv2(:,[184:258]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#DE7A22'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of conflict events','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

text(41.45,mean(ylim),'Oct. 29, 2018 to April 5, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)


subplot('Position',[0.565,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Ctv2(:,[184:258]),2);
MARo=sum(Ctv2(:,[109:183]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#DE7A22'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of conflict events','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;


text(41.45,mean(ylim),'April 6,2020 to Sept. 12, 2021','HorizontalAlignment','center','Rotation',90,'FontSize',20)

print(gcf,['Figure1_Conflict_COVID-19.png'],'-dpng','-r600');

%% Shelling

TestH=mean(Mt2(RC2==1,:),1);
TestG=mean(Mt2(RC2==0,:),1);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.047268907563025,0.624113475177305,0.942752100840336,0.365886524822696]); 

plot([109:258],TestH(109:258),'color',hex2rgb('#4C3F54'),'LineWidth',2); hold on; 
plot([109:258],TestG(109:258),'color',hex2rgb('#4C3F54'),'LineWidth',2,'LineStyle','-.'); hold on; 

plot(184.*ones(101,1),linspace(0,15,101),'k:','LineWidth',2);

startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(258-1)]],'mm/dd/yy');
xlim([109-0.5 258+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:258],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',18);
box off;
ylabel('Number of shellings and attacks','Fontsize',18);

text(-13.9,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');

legend('Houthi control','Government control')
legend boxoff;

subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Mt2(:,[109:183]),2);
MARo=sum(Mt2(:,[184:258]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#4C3F54'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of shellings and attacks','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

text(41.45,mean(ylim),'Oct. 29, 2018 to April 5, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)

subplot('Position',[0.565,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Mt2(:,[184:258]),2);
MARo=sum(Mt2(:,[109:183]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#4C3F54'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of shellings and attacks','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;


text(41.45,mean(ylim),'April 6,2020 to Sept. 12, 2021','HorizontalAlignment','center','Rotation',90,'FontSize',20)

print(gcf,['Figure1_Shelling-COVID-19.png'],'-dpng','-r600');

%% Diesel
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Dieselt2(RC2==1,:),1);
TestG=mean(Dieselt2(RC2==0,:),1);


plot([1:247],TestH,'color',[153,52,4]./255,'LineWidth',2); hold on; 
plot([1:247],TestG,'color',[153,52,4]./255,'LineWidth',2,'LineStyle','-.'); hold on; 


figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.047268907563025,0.624113475177305,0.942752100840336,0.365886524822696]); 

plot([120:247],TestH(120:247),'color',[153,52,4]./255,'LineWidth',2); hold on; 
plot([120:247],TestG(120:247),'color',[153,52,4]./255,'LineWidth',2,'LineStyle','-.'); hold on; 

plot(184.*ones(101,1),linspace(200,700,101),'k:','LineWidth',2);

startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(258-1)]],'mm/dd/yy');
xlim([120-0.5 247+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:258],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',18);
box off;
ylabel('Diesel price','Fontsize',18);


legend({'Houthi control','Government control'},'Location','NorthWest')
legend boxoff;

text(-13.9,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');


subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Dieselt2(:,[120:183]),2);
MARo=mean(Dieselt2(:,[184:247]),2);
MART=MAR(MAR>0);
MAR=(MAR-min(min(MARo(MARo>0)),min(MAR(MAR>0))))./(max(max(MARo),max(MAR))-min(min(MARo(MARo>0)),min(MAR(MAR>0))));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo(MARo>0)),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[153,52,4]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Average diesel price','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(41.45,mean(ylim),'Jan. 14, 2019 to April 5, 2020','HorizontalAlignment','center','Rotation',90,'FontSize',20)

subplot('Position',[0.565,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Dieselt2(:,[184:247]),2);
MARo=mean(Dieselt2(:,[120:183]),2);
MART=MAR(MAR>0);
MAR=(MAR-min(min(MARo(MARo>0)),min(MAR(MAR>0))))./(max(max(MARo),max(MAR))-min(min(MARo(MARo>0)),min(MAR(MAR>0))));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
load('Houthi_Hatching.mat','Stp','XtV');

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

for ii=1:length(S)
     if(RC(ii)==1)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            else
                mapshow(S(ii),'FaceColor','w','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',0); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end
     end
end


plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(min(MARo(MARo>0)),min(MART)),max(max(MARo),max(MART)),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[153,52,4]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Average diesel price','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

text(41.45,mean(ylim),'April 6, 2020 to June 27, 2021','HorizontalAlignment','center','Rotation',90,'FontSize',20)

print(gcf,['Figure1_Diesel-COVID-19.png'],'-dpng','-r600');