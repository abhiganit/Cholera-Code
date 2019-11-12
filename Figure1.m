close all;
clear;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,GNZI,maxtau] = LoadYemenData; % Load the data used to construct the figure

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

figure('units','normalized','outerposition',[0 0 1 1]);
load('Yemen_Gov_Incidence.mat')
IData=IData';
NW=153;
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
Pop=sum(PopS,1);
Test=10000.*sum(IData,1)./Pop;

subplot('Position',[0.05,0.6,0.945,0.39]); 
plot([1:153],Test,'color',hex2rgb('#C60000'),'Linewidth',2);
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
box off;
yh=ylabel({'Incidence per 10,000'},'Fontsize',18);
text(1.1*yh.Extent(1),0.99*max(ylim),'A','Fontsize',32,'FontWeight','bold');
subplot('Position',[0.05,0.16,0.945,0.39]); 
plot([1:153],sum(Rtv,1),'color',[5,112,176]./255,'Linewidth',2);
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
box off;
ylabel('Rainfall (mm)','Fontsize',18);

text(1.1*yh.Extent(1),0.99*max(ylim),'B','Fontsize',32,'FontWeight','bold');

figure('units','normalized','outerposition',[0 0 1 1]);

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
Pop=sum(PopS,1);
Test=10000.*sum(IData,1)./Pop;

subplot('Position',[0.05,0.6,0.945,0.39]); 
plot([1:153],sum(Ctv,1),'k','Linewidth',2); hold on;

plot([1:153],sum(Mt,1),'color',[0.65 0.65 0.65],'Linewidth',2);

legend({'Weekly conflict events','Air and shelling attacks'});
legend boxoff;
ylim([0 350]);
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',16,'Xminortick','on','YMinortick','on');
xtickangle(45);
box off;
ylabel('Number of events','Fontsize',18);
text(1.1*yh.Extent(1),0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.05,0.16,0.945,0.39]); 

plot([1:153],sum(tA,1),'color',[221,28,119]./255,'Linewidth',2);
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTick',[0:4]);
xtickangle(45);
box off;
ylim([0 4]);
ylabel('Attacks targeting water sources','Fontsize',18);

text(1.1*yh.Extent(1),0.99*max(ylim),'D','Fontsize',32,'FontWeight','bold');
xlabel('Week of report','Fontsize',18);



figure('units','normalized','outerposition',[0 0 1 1]);
%% Incidence
subplot('Position',[0.035,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
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
MAR=mean(WI,2);
MART=MAR;
MAR=MAR./max(MAR);
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(MART),100);
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
text(min(xlim),max(ylim)*0.98,'i)','Fontsize',24);

text(min(xlim)*0.98,max(ylim)*0.98,'E','Fontsize',32,'Fontweight','bold');
%% WASH
subplot('Position',[0.47,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(WPINm,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',[0 0.6 0.6],'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[0 0.6 0.6],'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Percentage of population in need of WaSH','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'ii)','Fontsize',24);
%% Food security
subplot('Position',[0.035,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(FPINm,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#2E4600'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Percentage of population in need of food security','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'iii)','Fontsize',24);
% %% Diesel

% 
% MAR=mean(Dieselt,2);
% MART=MAR(MAR>0);
% MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
% for ii=1:length(S)
%     if(MAR(ii)>=0)
%         mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
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
%     fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[153,52,4]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
%     if(rem(ii,5)==0)
%         h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);
%     end
% end
% text(mean(xlim),10.78,'Price of diesel','Fontsize',16,'HorizontalAlignment','center');
% ylim([11.7   19.0978]);
% 
% axis off;
% text(min(xlim),max(ylim)*0.98,'iv)','Fontsize',24);
%% Populatino Desnity
subplot('Position',[0.47,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(P,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',[136,86,167]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);

dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii),2)),'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[136,86,167]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'People per km^2','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'iv)','Fontsize',24);

figure('units','normalized','outerposition',[0 0 1 1]);

%% Hospital

subplot('Position',[0.035,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(H,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',[254,217,118]/255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[254,217,118]/255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Health facilities per 10,000','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'v)','Fontsize',24);
%% Targeted attacks

subplot('Position',[0.47,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=sum(tA,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',[221,28,119]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);
dA=min(MART):max(MART);
dX=linspace(min(xlim),max(xlim),length(dA));
ii=1;
h=text(dX(ii), 11.68,[num2str((dA(ii))) ],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:length(dA)
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[221,28,119]./255,'Facealpha',(ii-1)./length(dA),'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str((dA(ii)))],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of attacks on water systems','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'vi)','Fontsize',24);
%% Rainfall
subplot('Position',[0.035,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt


MAR=mean(Rtv,2);
MART=MAR(MAR>0);
MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
for ii=1:length(S)
    if(MAR(ii)>=0)
        mapshow(S(ii),'FaceColor',[5,112,176]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[5,112,176]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Rainfall (mm)','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'vii)','Fontsize',24);

% Legend for the location
subplot('Position',[0.47,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
for ii=1:length(S)
        mapshow(S(ii),'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2); hold on
        polyin = polyshape(S(ii).X,S(ii).Y);
        [x,y] = centroid(polyin);
        if(ii==19)
           text(44.453655324154,15.2224960733724, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
        elseif(ii==9)
            text(45.48567644535073,17.95860286252535, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
            annotation(gcf,'arrow',[0.590336134453782 0.558823529411765],[0.446821681864235 0.326241134751773]);
        elseif(ii==2)
             text(44.11595235285771,12.064063924702339, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
             annotation(gcf,'arrow',[0.558823529411765 0.573529411764706],[0.146909827760892 0.178318135764944]);
        else            
            text(x,y, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
        end
        if(ii<=8)
            text(45.45,12.8-0.3*(ii-1),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
        elseif(ii<=16)
            text(47.62,12.8-0.3*(ii-9),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
        else            
            text(50.505,12.8-0.3*(ii-17),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',12);
        end
end
box off;
xlim([41.7741   54.6472]);

ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'viii)','Fontsize',24);


figure('units','normalized','outerposition',[0 0 1 1]);
%% Conflict
subplot('Position',[0.035,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=sum(Ctv,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#FAAF08'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#FAAF08'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of conflict events','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'i)','Fontsize',24);

text(min(xlim)*0.98,max(ylim)*0.98,'E','Fontsize',32,'Fontweight','bold');
%% Shelling and air attacks
subplot('Position',[0.47,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=sum(Mt,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#4C3F54'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Total number of shelling and air attaks','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'ii)','Fontsize',24);
clear;