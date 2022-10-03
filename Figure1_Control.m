close all;
clear;
[WI2,~,tA2,Rtv2,~,P2,RC2,H2,WPINm2,FPINm2,Dieselt2,Wheatt2,V1,V2,GNZI,GV,maxtau,~,CI] = LoadYemenData; % Load the data used to construct the figure
RC=RC2;
WI2=WI2(GNZI,:);
tA2=tA2(GNZI,:);
Rtv2=Rtv2(GNZI,:);
P2=P2(GNZI,:);
RC2=RC2(GNZI);
H2=H2(GNZI);
WPINm2=WPINm2(GNZI,:);
FPINm2=FPINm2(GNZI,:);
Dieselt2=Dieselt2(GNZI,:);
Wheatt2=Wheatt2(GNZI,:);
% Need to increase the diesel price as we transformed it by subtracting the
% minimum
load('Diesel_Gov_Yemen.mat')
Dieselt2=Dieselt2+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt2=Wheatt2+min(Wheat(Wheat>0));

% Need to show the conflict without the transformation
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv2=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt2=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

Ctv2=Ctv2(GNZI,:);
Mt2=Mt2(GNZI,:);


S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


figure('units','normalized','outerposition',[0 0.1 1 0.6]);
%% Incidence
load('Yemen_Gov_Incidence.mat')
IData=IData';
IData=IData(GNZI,:);
NW=153;
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(GNZI,:);
PopH=sum(PopS(RC2==1,:),1);
PopG=sum(PopS(RC2==0,:),1);
TestH=10000.*sum(IData(RC2==1,:),1)./PopH;
TestG=10000.*sum(IData(RC2==0,:),1)./PopG;


subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 
plot([1:153],TestH,'color',hex2rgb('#C60000'),'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#C60000'),'MarkerEdgeColor',hex2rgb('#C60000'), 'MarkerIndices', 1:2:153); hold on;
plot([1:153],TestG,'color',hex2rgb('#C60000'),'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#C60000'), 'MarkerIndices', 1:2:153);
startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
yh=ylabel({'Incidence per 10,000'},'Fontsize',22);
legend({'Rebel control','Government control'},'Fontsize',20)
legend boxoff;
text(-21.107,0.99*max(ylim),'A','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

% Dot Hatching
% XS=linspace(41.7741,54.6472,71);
% YS=linspace(11.7,19.0978,71);
% [XSRCt,YSRCt]=meshgrid(XS,YS);
% XSRC=XSRCt(:);
% YSRC=YSRCt(:);
% 
% XS=[(XS(2:end)+XS(1:end-1))./2];
% YS=[(YS(2:end)+YS(1:end-1))./2];
% [XSRCt,YSRCt]=meshgrid(XS,YS);
% XSRC=[XSRC; XSRCt(:)];
% YSRC=[YSRC; YSRCt(:)];
% in=zeros(size(XSRC));
% for ii=1:length(S)
%     if(RC(ii)==1)
%         in=in+inpolygon(XSRC,YSRC,S(ii).X,S(ii).Y);
%     end
% end
% XSRC=XSRC(in>0);        
% YSRC=YSRC(in>0);




MAR=mean(WI2,2);
MART=MAR;
MAR=MAR./max(MAR);
for ii=1:length(S)
        if(ii~=21)
            if(ii<21)
                mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on    
            else
                mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on    
            end
        else
            mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
        end

%      if(RC(ii)==1)
%          if(ii==find(RC==1,1))
%              Stp=polyshape(S(ii).X,S(ii).Y);
%          else
%             Stp=union(Stp,polyshape(S(ii).X,S(ii).Y));
%          end
%      end
end

% Line Hatching
% 
% XS=linspace(41.7741,54.6472,1001);
% dXS=round(XS(2)-XS(1),4);
% YSb=linspace(11.7,19.0978,26);
% mx=0.5747;
% 
% XtV=[];
% for jj=1:26
%     in=inpolygon(XS,(YSb(jj)-mx.*XS(1))+mx.*XS,Stp.Vertices(:,1),Stp.Vertices(:,2));
%     XtS=XS(in);
%     LL=length(XtS);
%     xstart=1;
%     xend=1;
%     while(xend<LL)
%         while((round(XtS(xend+1)-XtS(xend),4)==dXS) && (xend<(length(XtS)-1)))
%             xend=xend+1;
%         end
%         plot(linspace(XtS(xstart),XtS(xend),101),(YSb(jj)-mx.*XS(1))+mx.*linspace(XtS(xstart),XtS(xend),101),'k:','LineWidth',1.1); hold on
%         XtV=[XtV; XtS(xstart) XtS(xend) (YSb(jj)-mx.*XS(1)) mx];
%         xstart=xend+1;
%         xend=xend+1;
%     end
%     
%     
%     in=inpolygon(XS,(YSb(jj)+mx.*XS(1))-mx.*XS,Stp.Vertices(:,1),Stp.Vertices(:,2));
%     XtS=XS(in);
%     LL=length(XtS);
%     xstart=1;
%     xend=1;
%     while(xend<LL)
%         while((round(XtS(xend+1)-XtS(xend),4)==dXS) && (xend<(length(XtS)-1)))
%             xend=xend+1;
%         end
%         plot(linspace(XtS(xstart),XtS(xend),101),(YSb(jj)+mx.*XS(1))-mx.*linspace(XtS(xstart),XtS(xend),101),'k:','LineWidth',1.1); hold on
%         XtV=[XtV; XtS(xstart) XtS(xend) (YSb(jj)+mx.*XS(1)) -mx];
%         xstart=xend+1;
%         xend=xend+1;
%     end
%     
% end
% save ('Rebel_Hatching.mat','Stp','XtV');
load('Rebel_Hatching.mat','Stp','XtV');

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


% scatter(XSRC,YSRC,5.25,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',20);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#C60000'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Incidence per 10,000','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);
axis off;
print(gcf,['Figure1_Incidence.png'],'-dpng','-r600');
%% Rainfall
TestH=mean(Rtv2(RC2==1,:),1);
TestG=mean(Rtv2(RC2==0,:),1);
figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

plot([1:153],TestH,'color',[5,112,176]./255,'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',[5,112,176]./255,'MarkerEdgeColor',[5,112,176]./255, 'MarkerIndices', 1:2:153); hold on;
plot([1:153],TestG,'color',[5,112,176]./255,'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',[5,112,176]./255, 'MarkerIndices', 1:2:153);

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
box off;
xlabel('Week of report','Fontsize',22);
ylabel('Rainfall (mm)','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20)
legend boxoff;
text(-21.107,40,'B','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Rtv2,2);
MART=MAR;
MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
for ii=1:length(S)
    if(ii~=21)
        if(ii<21)
            mapshow(S(ii),'FaceColor',[5,112,176]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[5,112,176]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',20);


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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[5,112,176]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Rainfall (mm)','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;

print(gcf,['Figure1_Rain.png'],'-dpng','-r600');

%% Conflict


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Ctv2(RC2==1,:),1);
TestG=mean(Ctv2(RC2==0,:),1);


figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

plot([1:153],TestH,'color',hex2rgb('#DE7A22'),'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#DE7A22'),'MarkerEdgeColor',hex2rgb('#DE7A22'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#DE7A22'),'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#DE7A22'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel('Number of conflict events','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(-21.107,0.99*max(ylim),'C','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=sum(Ctv2,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#DE7A22'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end

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


% scatter(XSRC,YSRC,5.25,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',20);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#DE7A22'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii))),'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Total number of conflict events','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Conflict.png'],'-dpng','-r600');
%% Shelling
figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

TestH=mean(Mt2(RC2==1,:),1);
TestG=mean(Mt2(RC2==0,:),1);

plot([1:153],TestH,'color',hex2rgb('#4C3F54'),'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#4C3F54'),'MarkerEdgeColor',hex2rgb('#4C3F54'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#4C3F54'),'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#4C3F54'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel('Number of shelling attacks','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(-21.107,15,'D','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt


MAR=sum(Mt2,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#4C3F54'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);


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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#4C3F54'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Total number of shelling and air attacks','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Shelling.png'],'-dpng','-r600');
%% Diesel
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
TestH=mean(Dieselt2(RC2==1,:),1);
TestG=mean(Dieselt2(RC2==0,:),1);


figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

plot([1:153],TestH,'color',[153,52,4]./255,'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',[153,52,4]./255,'MarkerEdgeColor',[153,52,4]./255, 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',[153,52,4]./255,'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',[153,52,4]./255, 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel('Diesel price','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(-21.107,0.99*max(ylim),'E','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Dieselt2,2);
MART=MAR;
MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);


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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[153,52,4]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Price of diesel','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Diesel.png'],'-dpng','-r600');
%% Wheat
figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

TestH=mean(Wheatt2(RC2==1,:),1);
TestG=mean(Wheatt2(RC2==0,:),1);

plot([1:153],TestH,'color',hex2rgb('#FAAF08'),'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#FAAF08'),'MarkerEdgeColor',hex2rgb('#FAAF08'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],TestG,'color',hex2rgb('#FAAF08'),'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#FAAF08'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel('Wheat price','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(-21.107,0.99*max(ylim),'F','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
MAR=mean(Wheatt2,2);
MART=MAR(MAR>0);
MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#FAAF08'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#FAAF08'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);

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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#FAAF08'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Price of wheat','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Wheat.png'],'-dpng','-r600');

%% target attacks


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019

figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

TestH=mean(tA2(RC2==1,:),1);
TestG=mean(tA2(RC2==0,:),1);
b=bar([1:153],[TestH; TestG]','stacked');

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel({'Number of attacks','on water systems'},'Fontsize',22);

ColorM=[[221,28,119]./255; [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=[221,28,119]./255;
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':'; 
legend('Rebel control','Government control');
legend boxoff;

text(-21.107,0.99*max(ylim),'G','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=sum(tA2,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[221,28,119]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[221,28,119]./255,'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=min(MART):max(MART);
dX=linspace(min(xlim),max(xlim),length(dA));
ii=1;
h=text(dX(ii), 11.68,[num2str((dA(ii))) ],'Rotation',270,'Fontsize',20);

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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:length(dA)
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[221,28,119]./255,'Facealpha',(ii-1)./length(dA),'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str((dA(ii)))],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Total number of attacks on water systems','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Attacks.png'],'-dpng','-r600');

%% WaSH
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP

figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 
load('WASH_PIN_Yemen.mat')
WPIN=WPIN(GNZI,:);
AP=AP(GNZI,:);
TestH=sum(WPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
TestG=sum(WPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
b=bar([2016 2017 2018 2019],[TestH;TestG]');

xlim([2015.5 2019.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',22,'YMinortick','on');
xtickangle(45);
xlabel('Year','Fontsize',22);
box off;
ylabel({'Proportion in need','of WaSH'},'Fontsize',22);

ColorM=[[0 0.6 0.6]; [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=[0 0.6 0.6];
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':';
ylim([0 0.75]);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(2014.92,0.99*max(ylim),'H','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(WPINm2,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',[0 0.6 0.6],'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',[0 0.6 0.6],'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',20);


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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[0 0.6 0.6],'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Percentage of population in need of WaSH','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_WASH.png'],'-dpng','-r600');

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
%% Food security
figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 
load('Food_PIN_Yemen.mat')
FPIN=FPIN(GNZI,:);
AP=AP(GNZI,:);
TestH=sum(FPIN(RC2==1,:),1)./sum(AP(RC2==1,:));
TestG=sum(FPIN(RC2==0,:),1)./sum(AP(RC2==0,:));
b=bar([2016 2017 2018 2019],[TestH;TestG]');

xlim([2015.5 2019.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[2016:2019],'Fontsize',22,'YMinortick','on');
xtickangle(45);
xlabel('Year','Fontsize',22);
box off;
ylabel({'Proportion in need','of food security'},'Fontsize',22);

ColorM=[hex2rgb('#2E4600'); [1 1 1]];
for ii=1:2
   b(ii).FaceColor=ColorM(ii,:); 
   b(ii).EdgeColor=hex2rgb('#2E4600');
   b(ii).LineWidth=2;
end

b(1).LineStyle='none'; 
b(2).LineStyle=':';


legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(2014.92,0.99*max(ylim),'I','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(FPINm2,2);
MART=MAR;
MAR=(MAR-min(MAR))./(max(MAR)-min(MAR));
for ii=1:length(S)
    if(ii~=21)        
        if(ii<21)
            mapshow(S(ii),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii)); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5,'FaceAlpha',MAR(ii-1)); hold on
        end
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',20);


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


% scatter(XSRC,YSRC,5.25,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#2E4600'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',20);
    end
end
text(mean(xlim),10.387196460176991,'Percentage of population in need of food security','Fontsize',22,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

plot(linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001),'k:','LineWidth',1.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Food.png'],'-dpng','-r600');
%% Vaccination
[V1,V2]= VaccinationTime(1,153);
V1=V1(GNZI,:); % Received at least one dose
TempH=cumsum(sum(V1(RC2==1,:),1))./PopH;
TempG=cumsum(sum(V1(RC2==0,:),1))./PopG;

figure('units','normalized','outerposition',[0 0.1 1 0.6]);
subplot('Position',[0.070903361344538,0.281879194630872,0.49062268907563,0.69221476510067]); 

plot([1:153],100.*TempH,'color',hex2rgb('#80BD9E'),'LineWidth',3,'LineStyle','-'); hold on; %,'Marker','s','MarkerFaceColor',hex2rgb('#80BD9E'),'MarkerEdgeColor',hex2rgb('#80BD9E'), 'MarkerIndices', 1:2:153); hold on; 
plot([1:153],100.*TempG,'color',hex2rgb('#80BD9E'),'LineWidth',2.25,'LineStyle',':'); %,'Marker','o','MarkerEdgeColor',hex2rgb('#80BD9E'), 'MarkerIndices', 1:2:153); hold on; 

startDateofSim = datenum('10-03-2016');% Start date
dW=8;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI2(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',20,'Xminortick','on','YMinortick','on');
xtickangle(45);
xlabel('Week of report','Fontsize',22);
box off;
ylabel('Vaccination uptake','Fontsize',22);

legend({'Rebel control','Government control'},'Fontsize',20,'Location','NorthWest')
legend boxoff;
text(-21.107,0.99*max(ylim),'J','Fontsize',32,'FontWeight','bold');

subplot('Position',[0.574054621848739,0.178378378378378,0.424663865546217,0.814414414414414]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

for ii=1:length(S)
        if(GV(ii)==0)
            mapshow(S(ii),'FaceColor','none','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5); hold on
        else
            mapshow(S(ii),'FaceColor',hex2rgb('#80BD9E'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5); hold on
        end
        polyin = polyshape(S(ii).X,S(ii).Y);
        [x,y] = centroid(polyin);
        if(ii==19)
           text(44.453655324154,15.2224960733724, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
        elseif(ii==9)
            text(45.48567644535073,17.95860286252535, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
            annotation(gcf,'arrow',[0.695903361344538 0.661764705882353],[0.836036036036036 0.609080204824885]);
        elseif(ii==2)
             text(44.11595235285771,12.064063924702339, num2str(ii),'HorizontalAlignment','center','Fontsize',12);
             annotation(gcf,'arrow',[0.658088235294118 0.676470588235293],[0.227027027027027 0.294723294723295]);
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

plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   


box off;
xlim([41.7741   54.6472]);

ylim([11.7   19.0978]);

plot(linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1),'k','LineWidth',2.5);
plot(linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1),'k','LineWidth',2.5);


plot(41.8.*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);
plot((41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001),'k','LineWidth',2.5);

text(41.7741+0.9,mean([18.4 18.95]),'Rebel Control','Fontsize',20);

axis off;
print(gcf,['Figure1_Vaccination.png'],'-dpng','-r600');