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



ph=0.21;
pw=0.573;
xx0=-0.1213;
xx1=0.2125;
xx2=0.548;

yy1=0.792;
yy2=0.53;
yy3=0.27;
yy4=0.0;
figure('units','normalized','outerposition',[0.1 0 0.5 0.8]);
subplot('Position',[xx0,yy1,pw,ph]);  % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

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
% save ('Houthi_Hatching.mat','Stp','XtV');
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


% scatter(XSRC,YSRC,5.25,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.675,num2str(dA(ii)),'Rotation',270,'Fontsize',8);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#C60000'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Incidence per 10,000','Fontsize',10,'HorizontalAlignment','center','Units','normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

text(0,0.93,'A','Fontsize',18,'FontWeight','bold','Units','Normalized');
axis off;

%% Rainfall
sbp=subplot('Position',[0.7,yy1,pw,ph]); 

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
h=text(dX(ii), 11.675,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',8);


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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(dA(ii),2))],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Rainfall (mm)','Fontsize',10,'HorizontalAlignment','center','Units','normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

temp_pos=sbp.Position;
sbp.Position=[xx1 temp_pos(2:end)];
text(0,0.93,'B','Fontsize',18,'FontWeight','bold','Units','Normalized');
axis off;

%% Conflict


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019

sbp=subplot('Position',[1.1,yy1,pw,ph]); 
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
h=text(dX(ii), 11.675,num2str(round(dA(ii))),'Rotation',270,'Fontsize',8);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],hex2rgb('#DE7A22'),'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,num2str(round(dA(ii))),'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Total number of conflict events','Fontsize',10,'HorizontalAlignment','center','Units','normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend
temp_pos=sbp.Position;
sbp.Position=[xx2 temp_pos(2:end)];
axis off;
text(0,0.93,'C','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Shelling


sbp=subplot('Position',[xx0,0.3,pw,ph]); 
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
h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);


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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Total number of shelling and air attacks','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

temp_pos=sbp.Position;
sbp.Position=[xx0 yy2 temp_pos(3:end)];

axis off;
text(0,0.93,'D','Fontsize',18,'FontWeight','bold','Units','Normalized');

%% Diesel
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019

sbp=subplot('Position',[0.7,0.3,pw,ph]); 

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
h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);


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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Price of diesel','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

temp_pos=sbp.Position;
sbp.Position=[xx1 yy2 temp_pos(3:end)];

text(0,0.93,'E','Fontsize',18,'FontWeight','bold','Units','Normalized');
axis off;
%% Wheat



sbp=subplot('Position',[1.1,0.3,pw,ph]); 
% Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
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
h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);

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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Price of wheat','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend
text(0,0.93,'F','Fontsize',18,'FontWeight','bold','Units','Normalized');


temp_pos=sbp.Position;
sbp.Position=[xx2 yy2 temp_pos(3:end)];

axis off;

%% target attacks


load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019


sbp=subplot('Position',[xx0,0.1,pw,ph]); 

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
h=text(dX(ii), 11.675,[num2str((dA(ii))) ],'Rotation',270,'Fontsize',8);

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
    if(rem(dA(ii),5)==0)
        h=text(dX(ii), 11.675,[num2str((dA(ii)))],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Total number of attacks on water systems','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend
text(0,0.93,'G','Fontsize',18,'FontWeight','bold','Units','Normalized');

temp_pos=sbp.Position;
sbp.Position=[xx0 yy3 temp_pos(3:end)];

axis off;

%% WaSH
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP

sbp=subplot('Position',[0.7 0.1,pw,ph]); 

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
h=text(dX(ii), 11.675,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',8);


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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Percentage of population in need of WaSH','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

text(0,0.93,'H','Fontsize',18,'FontWeight','bold','Units','Normalized');

temp_pos=sbp.Position;
sbp.Position=[xx1 yy3 temp_pos(3:end)];

axis off;

load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
%% Food security
 
load('Food_PIN_Yemen.mat')
FPIN=FPIN(GNZI,:);
AP=AP(GNZI,:);

sbp=subplot('Position',[1.1,0.1,pw,ph]); 
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
h=text(dX(ii), 11.675,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',8);


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
    if(rem(ii,10)==0)
        h=text(dX(ii), 11.675,[num2str(round(100.*dA(ii))) '%'],'Rotation',270,'Fontsize',8);
    end
end
text(0.5,-0.2,'Percentage of population in need of food security','Fontsize',10,'HorizontalAlignment','center','Units','Normalized');
ylim([11.7   19.0978]);

Houthi_Map_Legend

temp_pos=sbp.Position;
sbp.Position=[xx2 yy3 temp_pos(3:end)];

text(0,0.93,'I','Fontsize',18,'FontWeight','bold','Units','Normalized');
axis off;
%% Vaccination
sbp=subplot('Position',[xx0,-0.1,pw,ph]); 

[V1,V2]= VaccinationTime(1,153);
V1=V1(GNZI,:); % Received at least one dose

for ii=1:length(S)
    if(GV(ii)==0)
        mapshow(S(ii),'FaceColor','none','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5); hold on
    else
        mapshow(S(ii),'FaceColor',hex2rgb('#80BD9E'),'Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5); hold on
    end
end

for jj=1:length(XtV(:,1))
    plot(linspace(XtV(jj,1),XtV(jj,2),101),XtV(jj,3)+XtV(jj,4).*linspace(XtV(jj,1),XtV(jj,2),101),'k:','LineWidth',1.1); hold on
end

plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   


box off;
xlim([41.7741   54.6472]);

ylim([11.7   19.0978]);

Houthi_Map_Legend
text(0,0.93,'J','Fontsize',18,'FontWeight','bold','Units','Normalized');

temp_pos=sbp.Position;
sbp.Position=[xx0 yy4 temp_pos(3:end)];

axis off;

%% Map
sbp=subplot('Position',[0.7,-0.1,pw,ph]); 


for ii=1:length(S)
        mapshow(S(ii),'FaceColor','none','Edgecolor',[0.6 0.6 0.6],'LineWidth',1.5); hold on
        polyin = polyshape(S(ii).X,S(ii).Y);
        [x,y] = centroid(polyin);
        if(ii==19)
           text(44.453655324154,15.2224960733724, num2str(ii),'HorizontalAlignment','center','Fontsize',7);
        elseif(ii==9)
            text(45.48567644535073,17.95860286252535, num2str(ii),'HorizontalAlignment','center','Fontsize',7);
            annotation(gcf,'arrow',[0.432203389830508 0.407838983050847],[0.171206225680934 0.11284046692607]);
        elseif(ii==2)
             text(46.6732165503886,12.1097293568011, num2str(ii),'HorizontalAlignment','center','Fontsize',7);
             annotation(gcf,'arrow',[0.458686440677966 0.419491525423729],[0.012970168612192 0.0311284046692607]);
        else            
            text(x,y, num2str(ii),'HorizontalAlignment','center','Fontsize',7);
        end
        if(ii<=11)
            text(11+45.45,18.75-0.51*(ii-1),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',10);
        else
            text(15+47.62,18.75-0.51*(ii-12),[ num2str(ii) '. ' S(ii).ADM1_EN ],'Fontsize',10);
        end
end

plot(Stp,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2.5,'FaceAlpha',0); hold on   


box off;
xlim([41.7741   54.6472]);

ylim([11.7   19.0978]);

text(0,0.93,'K','Fontsize',18,'FontWeight','bold','Units','Normalized');

temp_pos=sbp.Position;
sbp.Position=[xx1 yy4 temp_pos(3:end)];

axis off;
print(gcf,['Figure2'],'-depsc','-r600');
print(gcf,['Figure2'],'-dpng','-r600');