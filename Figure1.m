close all;
clear;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,GNZI,maxtau] = LoadYemenData; % Load the data used to construct the figure

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

figure('units','normalized','outerposition',[0 0 1 1]);
%% Attack rate
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
    mapshow(S(ii),'FaceColor','r','Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],'r','Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,num2str(round(dA(ii),1)),'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Attack rate per 10,000','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'i)','Fontsize',16);
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
text(min(xlim),max(ylim)*0.98,'ii)','Fontsize',16);
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
text(min(xlim),max(ylim)*0.98,'iii)','Fontsize',16);
%% Diesel
subplot('Position',[0.47,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

MAR=mean(Dieselt,2);
MART=MAR(MAR>0);
MAR=(MAR-min(MAR(MAR>0)))./(max(MAR)-min(MAR(MAR>0)));
for ii=1:length(S)
    if(MAR(ii)>=0)
        mapshow(S(ii),'FaceColor',[153,52,4]./255,'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',MAR(ii)); hold on
    else
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.7 0.7 0.7],'LineWidth',2); hold on
    end
end
box off;
xlim([41.7741   54.6472]);
dA=linspace(min(MART),max(MART),100);
dX=linspace(min(xlim),max(xlim),100);
ii=1;
h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);

scatter(XSRC,YSRC,3,'k','filled');
for ii=2:100
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],[153,52,4]./255,'Facealpha',(ii-1)./99,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)))],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Price of diesel','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'iv)','Fontsize',16);
figure('units','normalized','outerposition',[0 0 1 1]);
%% Populatino Desnity
subplot('Position',[0.035,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

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
text(min(xlim),max(ylim)*0.98,'v)','Fontsize',16);
%% Hospital
subplot('Position',[0.47,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

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
text(min(xlim),max(ylim)*0.98,'vi)','Fontsize',16);
%% Targeted attacks
subplot('Position',[0.035,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

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
text(min(xlim),max(ylim)*0.98,'vii)','Fontsize',16);
%% Rainfall
subplot('Position',[0.47,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

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
text(mean(xlim),10.78,'Rainfall','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
text(min(xlim),max(ylim)*0.98,'viii)','Fontsize',16);
clear;