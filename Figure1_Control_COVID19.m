close all;
clear;
[WI2,~,tA2,Rtv2,~,P2,RC2,H2,WPINm2,FPINm2,~,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
WI2=WI2(GNZI,:);
tA2=tA2(GNZI,:);
Rtv2=Rtv2(GNZI,:);
P2=P2(GNZI,:);
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

[WI,~,tA,Rtv,~,P,RC,H,WPINm,FPINm,~,Wheatt,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
% Need to increase the diesel price as we transformed it by subtracting the
% minimum

Dieselt=DieselCOVID19;
load('Diesel_Gov_Yemen_COVID-19.mat');

Dieselt=Dieselt+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt=Wheatt+min(Wheat(Wheat>0));

% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

load('Conflict_COVID-19_Timeline.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,258);
load('Yemen_Air_Shelling_COVID-19.mat');
Mt=GLevelConflict(YASt,S,258); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

% figure('units','normalized','outerposition',[0 0 1 1]);
%% Incidence
% load('Yemen_Gov_Incidence.mat')
% IData=IData';
% IData=IData(GNZI,:);
% NW=153;
% load('PopulationSize_Yemen.mat');
% NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
% NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% % External effect due to IDP
% PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
% PopH=sum(PopS(RC2==1,:),1);
% PopG=sum(PopS(RC2==0,:),1);
% TestH=10000.*sum(IData(RC2==1,:),1)./PopH;
% TestG=10000.*sum(IData(RC2==0,:),1)./PopG;
% 
% Test=mean(TestH)./mean(TestG)
% 
% 
% Test=TestH./TestG;
% 
% Test=Test(Test~=Inf);
% 
% Test=Test(~isnan(Test));
% mean(Test)
% median(Test)
% 
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% plot([1:153],TestH,'color',hex2rgb('#C60000'),'LineWidth',2); hold on;
% plot([1:153],TestG,'color',hex2rgb('#C60000'),'LineWidth',2,'LineStyle','-.');
% startDateofSim = datenum('10-03-2016');% Start date
% dW=5;
% XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% xlim([0.5 length(WI2(1,:))+0.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
% xtickangle(45);
% xlabel('Week of report','Fontsize',18);
% box off;
% yh=ylabel({'Incidence per 10,000'},'Fontsize',18);
% legend('Rebel control','Government control')
% legend boxoff;
% text(-13.9,0.99*max(ylim),'A','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
XS=linspace(41.7741,54.6472,101);
YS=linspace(11.7,19.0978,101);
[XSRCt,YSRCt]=meshgrid(XS,YS);
XSRC=XSRCt(:);
YSRC=YSRCt(:);

XS=[(XS(2:end)+XS(1:end-1))./2];
YS=[(YS(2:end)+YS(1:end-1))./2];
[XSRCt,YSRCt]=meshgrid(XS,YS);
XSRC=[XSRC; XSRCt(:)];
YSRC=[YSRC; YSRCt(:)];
in=zeros(size(XSRC));
for ii=1:length(S)
    if(RC(ii)==1)
        in=in+inpolygon(XSRC,YSRC,S(ii).X,S(ii).Y);
    end
end
XSRC=XSRC(in>0);        
YSRC=YSRC(in>0);


figure('units','normalized','outerposition',[0 0 1 1]);
load('PopulationSize_Yemen.mat');
load('CholeraIncidence_COVID-19.mat');
CC19=cell2mat(CholeraCOVID19(:,2:end))';
CC19Date=datenum(CholeraCOVID19(:,1));
% CC19=CC19(GNZI,:);

% AP=AP(GNZI,:);
PopH=sum(AP(RC==1,4));
PopG=sum(AP(RC==0,4));
TestH=10000.*sum(CC19(RC==1,:),1)./PopH;
TestG=10000.*sum(CC19(RC==0,:),1)./PopG;


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
legend('Rebel control','Government control')
legend boxoff;

subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=10000.*mean(CC19(:,1:55)./repmat(AP(:,4),1,55),2);
MARo=10000.*mean(CC19(:,56:end)./repmat(AP(:,4),1,length(CC19(1,56:end))),2);
MART=MAR;
MAR=MAR./max(max(MAR),max(MARo));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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
    mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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
% %% Rainfall
% TestH=mean(Rtv2(RC2==1,:),1);
% TestG=mean(Rtv2(RC2==0,:),1);
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% 
% plot([1:153],TestH,'color',[5,112,176]./255,'LineWidth',2); hold on;
% plot([1:153],TestG,'color',[5,112,176]./255,'LineWidth',2,'LineStyle','-.');
% 
% startDateofSim = datenum('10-03-2016');% Start date
% dW=5;
% XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% xlim([0.5 length(WI2(1,:))+0.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
% xtickangle(45);
% box off;
% xlabel('Week of report','Fontsize',18);
% ylabel('Rainfall (mm)','Fontsize',18);
% 
% text(-13.9,0.99*max(ylim),'B','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% MAR=mean(Rtv,2);
% MART=MAR(MAR>0);
% MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
% for ii=1:length(S)
%     if(MAR(ii)>=0)
%         mapshow(S(ii),'FaceColor',[5,112,176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
%     else
%         mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
%     end
% end
% box off;
% xlim([41.7741   54.6472]);
% dA=linspace(min(MART),max(MART),100);
% dX=linspace(min(xlim),max(xlim),100);
% ii=1;
% h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);
% 
% scatter(XSRC,YSRC,3,'k','filled');
% for ii=2:100
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[5,112,176]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Rainfall (mm)','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% 
% print(gcf,['Figure1_Rain.png'],'-dpng','-r600');

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
ylabel('Number of conflict evetns','Fontsize',18);

text(-13.9,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');


subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Ctv(:,[109:183]),2);
MARo=sum(Ctv(:,[184:258]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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
MAR=sum(Ctv(:,[184:258]),2);
MARo=sum(Ctv(:,[109:183]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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


subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Mt(:,[109:183]),2);
MARo=sum(Mt(:,[184:258]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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
MAR=sum(Mt(:,[184:258]),2);
MARo=sum(Mt(:,[109:183]),2);
MART=MAR;
MAR=(MAR-min(min(MARo),min(MAR)))./(max(max(MARo),max(MAR))-min(min(MARo),min(MAR)));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
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

text(-13.9,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');


subplot('Position',[0.017268907563025,0.08,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Dieselt(:,[120:183]),2);
MARo=mean(Dieselt(:,[184:247]),2);
MART=MAR(MAR>0);
MAR=(MAR-min(min(MARo(MARo>0)),min(MAR(MAR>0))))./(max(max(MARo),max(MAR))-min(min(MARo(MARo>0)),min(MAR(MAR>0))));
for ii=1:length(S)
    if(MAR(ii)>=0)
        mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
scatter(XSRC,YSRC,3,'k','filled');
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
MAR=mean(Dieselt(:,[184:247]),2);
MARo=mean(Dieselt(:,[120:183]),2);
MART=MAR(MAR>0);
MAR=(MAR-min(min(MARo(MARo>0)),min(MAR(MAR>0))))./(max(max(MARo),max(MAR))-min(min(MARo(MARo>0)),min(MAR(MAR>0))));
for ii=1:length(S)
    if(MAR(ii)>=0)
        mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
scatter(XSRC,YSRC,3,'k','filled');
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
% %% Wheat
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% 
% TestH=mean(Wheatt2(RC2==1,:),1);
% TestG=mean(Wheatt2(RC2==0,:),1);
% 
% plot([1:153],TestH,'color',hex2rgb('#FAAF08'),'LineWidth',2); hold on; 
% plot([1:153],TestG,'color',hex2rgb('#FAAF08'),'LineWidth',2,'LineStyle','-.'); hold on; 
% 
% startDateofSim = datenum('10-03-2016');% Start date
% dW=5;
% XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% xlim([0.5 length(WI2(1,:))+0.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
% xtickangle(45);
% xlabel('Week of report','Fontsize',18);
% box off;
% ylabel('Wheat price','Fontsize',18);
% 
% text(-13.9,0.99*max(ylim),'F','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% MAR=mean(Wheatt,2);
% MART=MAR(MAR>0);
% MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
% for ii=1:length(S)
%     if(MAR(ii)>=0)
%         mapshow(S(ii),'FaceColor',hex2rgb('#FAAF08'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
%     else
%         mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
%     end
% end
% box off;
% xlim([41.7741   54.6472]);
% dA=linspace(min(MART),max(MART),100);
% dX=linspace(min(xlim),max(xlim),100);
% ii=1;
% h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);
% 
% scatter(XSRC,YSRC,3,'k','filled');
% for ii=2:100
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#FAAF08'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Price of wheat','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% print(gcf,['Figure1_Wheat.png'],'-dpng','-r600');
% 
% %% target attacks
% 
% 
% load('PopulationSize_Yemen.mat');
% NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
% NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% 
% TestH=mean(tA2(RC2==1,:),1);
% TestG=mean(tA2(RC2==0,:),1);
% b=bar([1:153],[TestH; TestG]','stacked');
% 
% startDateofSim = datenum('10-03-2016');% Start date
% dW=5;
% XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% xlim([0.5 length(WI2(1,:))+0.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
% xtickangle(45);
% xlabel('Week of report','Fontsize',18);
% box off;
% ylabel('Number of attacks','Fontsize',18);
% 
% ColorM=[[221,28,119]./255; [1 1 1]];
% for ii=1:2
%    b(ii).FaceColor=ColorM(ii,:); 
%    b(ii).EdgeColor=[221,28,119]./255;
%    b(ii).LineWidth=2;
% end
% 
% b(1).LineStyle='none'; 
% b(2).LineStyle='-.'; 
% legend('Rebel control','Government control');
% legend boxoff;
% 
% text(-13.9,0.99*max(ylim),'G','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% 
% MAR=sum(tA,2);
% MART=MAR;
% MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
% for ii=1:length(S)
%     mapshow(S(ii),'FaceColor',[221,28,119]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
% end
% box off;
% xlim([41.7741   54.6472]);
% dA=min(MART):max(MART);
% dX=linspace(min(xlim),max(xlim),length(dA));
% ii=1;
% h=text(dX(ii), 11.68,[num2str((dA(ii))) ],'Rotation',270,'Fontsize',14);
% 
% scatter(XSRC,YSRC,3,'k','filled');
% for ii=2:length(dA)
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[221,28,119]./255,'Facealpha',(ii-1)./length(dA),'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str((dA(ii)))],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Total number of attacks on water systems','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% print(gcf,['Figure1_Attacks.png'],'-dpng','-r600');
% 
% %% WaSH
% load('PopulationSize_Yemen.mat');
% NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
% NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% % External effect due to IDP
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% load('WASH_PIN_Yemen.mat')
% WPIN=WPIN(GNZI,:);
% AP=AP(GNZI,:);
% TestH=sum(WPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
% TestG=sum(WPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
% b=bar([2016 2017 2018 2019],[TestH;TestG]');
% 
% xlim([2015.5 2019.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',16,'YMinortick','on');
% xtickangle(45);
% xlabel('Year','Fontsize',18);
% box off;
% ylabel('Proportion in need of WaSH','Fontsize',18);
% 
% ColorM=[[0 0.6 0.6]; [1 1 1]];
% for ii=1:2
%    b(ii).FaceColor=ColorM(ii,:); 
%    b(ii).EdgeColor=[0 0.6 0.6];
%    b(ii).LineWidth=2;
% end
% 
% b(1).LineStyle='none'; 
% b(2).LineStyle='-.';
% 
% text(2015.14,0.99*max(ylim),'H','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% 
% MAR=mean(WPINm,2);
% MART=MAR;
% MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
% for ii=1:length(S)
%     mapshow(S(ii),'FaceColor',[0 0.6 0.6],'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
% end
% box off;
% xlim([41.7741   54.6472]);
% dA=linspace(min(MART),max(MART),100);
% dX=linspace(min(xlim),max(xlim),100);
% ii=1;
% h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
% 
% scatter(XSRC,YSRC,3,'k','filled');
% for ii=2:100
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[0 0.6 0.6],'Facealpha',(ii-1)./99,'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Percentage of population in need of WaSH','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% print(gcf,['Figure1_WASH.png'],'-dpng','-r600');
% load('PopulationSize_Yemen.mat');
% NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
% NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% %% Food security
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% load('Food_PIN_Yemen.mat')
% FPIN=FPIN(GNZI,:);
% AP=AP(GNZI,:);
% TestH=sum(FPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
% TestG=sum(FPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
% b=bar([2016 2017 2018 2019],[TestH;TestG]');
% 
% xlim([2015.5 2019.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',16,'YMinortick','on');
% xtickangle(45);
% xlabel('Year','Fontsize',18);
% box off;
% ylabel('Proportion in need of food security','Fontsize',18);
% 
% ColorM=[hex2rgb('#2E4600'); [1 1 1]];
% for ii=1:2
%    b(ii).FaceColor=ColorM(ii,:); 
%    b(ii).EdgeColor=hex2rgb('#2E4600');
%    b(ii).LineWidth=2;
% end
% 
% b(1).LineStyle='none'; 
% b(2).LineStyle='-.';
% 
% text(2015.15,0.99*max(ylim),'I','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% 
% MAR=mean(FPINm,2);
% MART=MAR;
% MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
% for ii=1:length(S)
%     mapshow(S(ii),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
% end
% box off;
% xlim([41.7741   54.6472]);
% dA=linspace(min(MART),max(MART),100);
% dX=linspace(min(xlim),max(xlim),100);
% ii=1;
% h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
% 
% scatter(XSRC,YSRC,3,'k','filled');
% for ii=2:100
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#2E4600'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Percentage of population in need of food security','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% print(gcf,['Figure1_Food.png'],'-dpng','-r600');
% %% Vaccination
% [V1,V2]= VaccinationTime(1,153);
% V1=V1(GNZI,:); % Received at least one dose
% TempH=cumsum(sum(V1(RC2==1,:),1))./PopH;
% TempG=cumsum(sum(V1(RC2==0,:),1))./PopG;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot('Position',[0.047268907563025,0.603900000000001,0.511631092436975,0.3861]); 
% 
% plot([1:153],100.*TempH,'color',hex2rgb('#80BD9E'),'LineWidth',2); hold on; 
% plot([1:153],100.*TempG,'color',hex2rgb('#80BD9E'),'LineWidth',2,'LineStyle','-.'); hold on; 
% 
% startDateofSim = datenum('10-03-2016');% Start date
% dW=5;
% XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% xlim([0.5 length(WI2(1,:))+0.5]);
% set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YMinortick','on');
% xtickangle(45);
% xlabel('Week of report','Fontsize',18);
% box off;
% ylabel('Vaccination uptake','Fontsize',18);
% 
% text(-13.9,0.99*max(ylim),'J','Fontsize',32,'FontWeight','bold');
% 
% subplot('Position',[0.565,0.55,0.45,0.45]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
% 
% for ii=1:length(S)
%         if(GV(ii)==0)
%             mapshow(S(ii),'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2); hold on
%         else
%             mapshow(S(ii),'FaceColor',hex2rgb('#80BD9E'),'Edgecolor',[0 0 0],'LineWidth',2); hold on
%         end
%         polyin = polyshape(S(ii).X,S(ii).Y);
%         [x,y] = centroid(polyin);
%         if(ii==19)
%            text(44.453655324154,15.2224960733724, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
%         elseif(ii==9)
%             text(45.48567644535073,17.95860286252535, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
%             annotation(gcf,'arrow',[0.701680672268907 0.668067226890756],[0.918959473150965 0.789260385005066]);
%         elseif(ii==2)
%              text(44.11595235285771,12.064063924702339, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
%              annotation(gcf,'arrow',[0.664390756302521 0.682247899159664],[0.578520770010134 0.619047619047619]);
%         else            
%             text(x,y, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
%         end
%         if(ii<=8)
%             text(45.45,12.8-0.3*(ii-1),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
%         elseif(ii<=16)
%             text(47.62,12.8-0.3*(ii-9),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
%         else            
%             text(50.505,12.8-0.3*(ii-17),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
%         end
% end
% box off;
% xlim([41.7741   54.6472]);
% 
% ylim([11.7   19.0978]);
% 
% axis off;
% print(gcf,['Figure1_Vaccination.png'],'-dpng','-r600');