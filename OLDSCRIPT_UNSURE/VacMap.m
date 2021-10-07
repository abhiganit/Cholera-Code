close all;
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
[~,~,~,~,~,~,RC,~,~,~,~,~,VT1,VT2,~,~,~] = LoadYemenData;
VT=sum(VT1,2);
FA=VT./max(VT);

figure('units','normalized','outerposition',[0 0 1 1]);
%% At least one dose
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
for ii=1:length(S)
    mapshow(S(ii),'FaceColor','m','Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',FA(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,100.*max(VT),40);
dX=linspace(min(xlim),max(xlim),40);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:40
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],'m','Facealpha',(ii-1)./39,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii))) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Coverage for at least one dose of OCV','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

%% At least two doses
VT=sum(VT2,2);
FA=VT./max(VT);
subplot('Position',[0.47,0.6,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

for ii=1:length(S)
    mapshow(S(ii),'FaceColor','c','Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',FA(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,100.*max(VT),40);
dX=linspace(min(xlim),max(xlim),40);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:40
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],'c','Facealpha',(ii-1)./39,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii))) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Coverage for at least two doses of OCV','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;
