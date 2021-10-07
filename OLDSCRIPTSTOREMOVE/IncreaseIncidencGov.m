% clc;
% clear;
LW=8;
WeekInt=90;
% NS=10^3;
% load('Yemen_Gov_Incidence.mat')
% IData=IData(:,[1:20 22])';
% RA=zeros(21,NS);
% RB=zeros(21,NS);
% IDataT=IData(:,(WeekInt):(WeekInt+(LW-1)));
% IDataTB=IData(:,(WeekInt-1-(LW-1)):(WeekInt-1));
% for ii=1:NS
%     if(ii==1)
%         IDataS=IDataT;
%         IDataSB=IDataTB;
%     else
%         IDataS=poissrnd(IDataT);
%         IDataSB=poissrnd(IDataTB);
%     end
%     IDataS=log(10^(-8)+IDataS);
%     IDataSB=log(10^(-8)+IDataSB);
%     for jj=1:21
%        temp=fmincon(@(x)sum(((x(1).*[0:(LW-1)]+x(2))-IDataS(jj,:)).^2),[1,0.1]);
%        RA(jj,ii)=temp(1);
%        temp=fmincon(@(x)sum(((x(1).*[0:(LW-1)]+x(2))-IDataSB(jj,:)).^2),[1,0.1]);
%        RB(jj,ii)=temp(1);
%     end
% end
% clc
% h2=zeros(21,1);
% p2=zeros(21,1);
% for ii=1:21
%   [h2(ii),p2(ii)]=kstest2(RA(ii,:),RB(ii,:));
% end
% d=mean(RA,2)-mean(RB,2);
% h2(isnan(h2))=0;
% sum(h2);
close all;
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

XGL={S([1:20 22]).ADM1_EN};
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.066176470588235,0.282722843934367,0.929621848739496,0.700027915944055]);
VP=[RB(1,:)' RA(1,:)'];
for mm=2:21
   VP=[VP RB(mm,:)' RA(mm,:)'];
end
v=violinplot(VP); hold on;
 for ii=1:21
    v(1+2.*(ii-1)).ViolinColor='b'; 
    v(2.*ii).ViolinColor='r'; 
%     v(1+2.*(ii-1)).MedianColor='b';
%     v(2.*ii).MedianColor='r';
    v(1+2.*(ii-1)).BoxColor='b';
    v(2.*ii).BoxColor='r';
    v(1+2.*(ii-1)).EdgeColor='b';
    v(2.*ii).EdgeColor='r';
    plot([1+2.*(ii-1) 2*ii], [median(VP(:,1+2.*(ii-1))) median(VP(:,2*ii))],'k','LineWidth',2);
    if(p2(ii)<0.05)
        if(d(ii)>0)
            scatter(mean([1+2.*(ii-1) 2*ii]), 4.6,40,[hex2rgb('#C60000')],'filled');
        else            
            scatter(mean([1+2.*(ii-1) 2*ii]), 4.6,40,[0.7 0.7 0.7],'filled');
        end
    end
 end
 xlim([0 43]);
 set(gca,'LineWidth',2,'tickdir','out','Xtick',[1.5:2:41.5],'YTick',[-4.5:0.5:4.5],'Fontsize',18,'Xticklabel',XGL)
 xtickangle(45)
 xlabel('Governorate','Fontsize',22)
 ylim([-4.5 4.7]);
 yh=ylabel('Slope of log incidence','Fontsize',22);
 text(yh.Extent(1),max(ylim),'A','Fontsize',32,'FontWeight','bold');
 text(3,-4,'B','Fontsize',32,'FontWeight','bold');
 text(5,-4,'C','Fontsize',32,'FontWeight','bold');
 startDateofSim = datenum('10-03-2016');% Start date
 astart=datestr(startDateofSim+7.*(WeekInt-1)); % subtract one from the week of interest because of the indexin (i.e. if one we want the start data);
 aend=datestr(startDateofSim+7.*(WeekInt+(LW)-1)-1); % subtract one from the week of interest because of the indexin (i.e. if one we want the start data); do not subtract one from lw as we weill be producing the end of the week. that is why we subtratc one outside the multiplication by 7
 
 bstart=datestr(startDateofSim+7.*((WeekInt-1-(LW-1))-1)); % subtract one from the week of interest because of the indexin (i.e. if one we want the start data);
 bend=datestr(startDateofSim+7.*(WeekInt-1)-1); % subtract one from the week of interest because of the indexin (i.e. if one we want the start data); do not subtract one from lw as we weill be producing the end of the week. that is why we subtratc one outside the multiplication by 7
 text(36.30960451977403,4.377,[astart ' to ' aend],'Color','r','Fontsize',16);
 text(36.30960451977403,4.06,[bstart ' to ' bend],'Color','b','Fontsize',16);
 
 
 [~,~,~,~,~,~,RC,~,~,~,~,~,~,~,~,~,~,~,~] = LoadYemenData;
 S=S([1:20 22]);
 RC=RC([1:20 22]);
figure('units','normalized','outerposition',[0 0 1 1]);
 subplot('Position',[0.02,0.02,0.98,0.98]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
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
    if(d(ii)*h2(ii)>0)
        mapshow(S(ii),'FaceColor',hex2rgb('#C60000'),'Edgecolor',[0 0 0],'LineWidth',2); hold on 
    elseif (d(ii)*h2(ii)<0)
        mapshow(S(ii),'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0 0 0],'LineWidth',2); hold on 
    else
        mapshow(S(ii),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',2); hold on 
    end
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);ylim([11.7   19.0978]);

axis off;
text(43.090731334022756,18.648074663909004,'Siginifcant increase in incidence rate','Color',hex2rgb('#C60000'),'Fontsize',22);
text(43.090731334022756,18.426,'Siginifcant decrease in incidence rate','Color',[0.7 0.7 0.7],'Fontsize',22);

figure('units','normalized','outerposition',[0 0 1 1]);
 subplot('Position',[0.02,0.02,0.98,0.98]); 
HPD=[3 5 6 8 9 10 11 13 14 17 18 19 21];
IND=zeros(21,1);
IND(HPD)=1;
for ii=1:length(S)
    if(IND(ii)>0)
        mapshow(S(ii),'FaceColor',hex2rgb('#DFE166'),'Edgecolor',[0 0 0],'LineWidth',2); hold on 
    else
        mapshow(S(ii),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',2); hold on 
    end
end
scatter(XSRC,YSRC,3,'k','filled');
box off;
xlim([41.7741   54.6472]);ylim([11.7   19.0978]);

axis off;
text(43.090731334022756,18.648074663909004,'Dependency on ports in Al Hodeidah','Color',hex2rgb('#DFE166'),'Fontsize',22);

