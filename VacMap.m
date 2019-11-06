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

%% Determine the maximum effectiveness

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
PDS=0.90;
atest=0;
%% Forward selection
load(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);


% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF,mln,a,KV,dV]=RetParameterPS(parv(end,:),XUv(end,:));

dV1=ImpactAttack(V1,0,dV,2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2,0,dV,2,maxtau);  % Two week delay until acquire immunity

EOVC=EffectOCV(dV1,KV(1),dV2,0.*KV(2));

MVE=max(EOVC')';

FA=MVE./max(MVE);
subplot('Position',[0.035,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

for ii=1:length(S)
    mapshow(S(ii),'FaceColor','c','Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',FA(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,100.*max(MVE),40);
dX=linspace(min(xlim),max(xlim),40);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:40
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],'c','Facealpha',(ii-1)./39,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii))) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Maximim effectiveness for at least one dose','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;

subplot('Position',[0.47,0.12,0.4,0.4]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

EOVC=EffectOCV(dV1,KV(1),dV2,KV(2));

MVE=max(EOVC')';

FA=MVE./max(MVE);

for ii=1:length(S)
    mapshow(S(ii),'FaceColor','c','Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',FA(ii)); hold on    
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);
dA=linspace(0,100.*max(MVE),40);
dX=linspace(min(xlim),max(xlim),40);
ii=1;
h=text(dX(ii), 11.68,num2str(dA(ii)),'Rotation',270,'Fontsize',14);

for ii=2:40
    fill([dX(ii-1) dX(ii-1) dX(ii) dX(ii)],[11.7 12 12 11.7],'c','Facealpha',(ii-1)./39,'Edgealpha',0);    
    if(rem(ii,5)==0)
        h=text(dX(ii), 11.68,[num2str(round(dA(ii)),6) '%'],'Rotation',270,'Fontsize',14);
    end
end
text(mean(xlim),10.78,'Overall effectiveness','Fontsize',16,'HorizontalAlignment','center');
ylim([11.7   19.0978]);

axis off;