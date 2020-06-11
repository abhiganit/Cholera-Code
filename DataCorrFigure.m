close all;


[~,~,tA,~,~,~,~,~,WPINm,FPINm,Dieselt,Wheatt,~,~,GNZI,~,~,~,~] = LoadYemenData;

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

% Need to increase the diesel price as we transformed it by subtracting the
% minimum
load('Diesel_Gov_Yemen.mat')
Dieselt=Dieselt+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
Wheatt=Wheatt+min(Wheat(Wheat>0));

D=Dieselt(GNZI,:);
W=Wheatt(GNZI,:);
C=Ctv(GNZI,:);
S=Mt(GNZI,:);
Ws=WPINm(GNZI,:);
Fs=FPINm(GNZI,:);
TA=tA(GNZI,:);

WsgU=zeros(21,4);
FsgU=zeros(21,4);
CgU=zeros(21,4);
for ii=1:21
   WsgU(ii,:)=unique(Ws(ii,:)); 
   for jj=1:4
    f=find(Ws(ii,:)==WsgU(ii,jj));
    CgU(ii,jj)=mean(C(ii,f));
   end
end

[r,p]=corr(CgU(:),WsgU(:));

for ii=1:21
   FsgU(ii,:)=unique(Fs(ii,:)); 
   for jj=1:4
    f=find(Fs(ii,:)==FsgU(ii,jj));
    CgU(ii,jj)=mean(C(ii,f));
   end
end

[r,p]=corr(CgU(:),FsgU(:));

C=sum(C,1)';
S=sum(S,1)';
D=mean(D,1)';
W=mean(W,1)';
Ws=mean(Ws,1)';
Fs=mean(Fs,1)';
TA=sum(TA,1)';

X=[C S TA D W Ws Fs]';

r=zeros(7,7);

for xx=1:7
    for yy=1:7
        rt=zeros(70,1);
        pt=zeros(70,1);
        for ii=1:70
            T=X(yy,ii:end);
            for jj=(ii-1):-1:1
                T=T+X(yy,jj:end-(ii-jj));
            end
            [rt(ii),pt(ii)]=corr(T',X(xx,ii:end)');
        end
        r(xx,yy)=max([rt(pt<0.05);rt(pt==min(pt))]);
    end 
end
rs=r;
dx=[1:6];
rs(rs<0)=1;
for ii=2:6
    f=find(r>=0.2*(ii-2));
    g=find(r(f)<0.2*(ii-1));
    f=f(g);
    rs(f)=ii;
end
rs(r==1)=6;
imagesc(rs)
box off;
set(gca,'Tickdir','out','LineWidth',2,'XTickLabel',{'Conflict','Sheliing','Targeted','Diesel','Wheat','WaSH','Malnutrition'},'YTickLabel',{'Conflict','Sheliing','Targeted','Diesel','Wheat','WaSH','Malnutrition'});
xtickangle(45);
xlabel('Cumulative effect');
ylabel('Weekly value');
h=colorbar;
x=[ones(1,3);hex2rgb('#C4DFE6');hex2rgb('#66A5Ad');hex2rgb('#07575B');hex2rgb('#003B46')];
colormap([ hex2rgb('#D5D6D2');x]); 
caxis([0.5 6.5]);
h.YTick=[0.5:1:6.5];
h.YTickLabel={'-1','0','0.2','0.4','0.6','0.8','1'};
h.TickDirection='out';
% %% Diesel conflict
% r=zeros(71,1);
% for ii=1:71
% T=C(ii:end);
% for jj=(ii-1):-1:1
% T=T+C(jj:end-(ii-jj));
% end
% r(ii)=corr(T,D(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of conflict','Fontsize',22)
% box off;
% title('Diesel and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(C(1:(end-(ii-1))),D(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of conflict','Fontsize',22)
% box off;
% title('Diesel and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% %% Diesel and shellings
% 
% r=zeros(71,1);
% for ii=1:71
% T=S(ii:end);
% for jj=(ii-1):-1:1
% T=T+S(jj:end-(ii-jj));
% end
% r(ii)=corr(T,D(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shellings','Fontsize',22)
% box off;
% title('Diesel and Shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(S(1:(end-(ii-1))),D(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of shellings','Fontsize',22)
% box off;
% title('Diesel and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% wheat and conflict
% r=zeros(71,1);
% for ii=1:71
% T=C(ii:end);
% for jj=(ii-1):-1:1
% T=T+C(jj:end-(ii-jj));
% end
% r(ii)=corr(T,W(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of conflict','Fontsize',22)
% box off;
% title('Wheat and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(C(1:(end-(ii-1))),W(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of conflict','Fontsize',22)
% box off;
% title('Wheat and conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% %% wheat and shellings
% 
% r=zeros(71,1);
% for ii=1:71
% T=S(ii:end);
% for jj=(ii-1):-1:1
% T=T+S(jj:end-(ii-jj));
% end
% r(ii)=corr(T,W(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shellings','Fontsize',22)
% box off;
% title('Wheat and Shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(S(1:(end-(ii-1))),W(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of shellings','Fontsize',22)
% box off;
% title('Wheat and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% WaSH conflict
% r=zeros(71,1);
% for ii=1:71
% T=C(ii:end);
% for jj=(ii-1):-1:1
% T=T+C(jj:end-(ii-jj));
% end
% r(ii)=corr(T,Ws(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of conflict','Fontsize',22)
% box off;
% title('PIN WaSH and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(C(1:(end-(ii-1))),Ws(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of conflict','Fontsize',22)
% box off;
% title('PIN WaSH and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% WaSH shellings
% r=zeros(71,1);
% for ii=1:71
% T=S(ii:end);
% for jj=(ii-1):-1:1
% T=T+S(jj:end-(ii-jj));
% end
% r(ii)=corr(T,Ws(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shellings','Fontsize',22)
% box off;
% title('PIN WaSH and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(S(1:(end-(ii-1))),Ws(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of shellings','Fontsize',22)
% box off;
% title('PIN WaSH and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% Food conflict
% r=zeros(71,1);
% for ii=1:71
% T=C(ii:end);
% for jj=(ii-1):-1:1
% T=T+C(jj:end-(ii-jj));
% end
% r(ii)=corr(T,Fs(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of conflict','Fontsize',22)
% box off;
% title('PIN Food and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(C(1:(end-(ii-1))),Fs(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of conflict','Fontsize',22)
% box off;
% title('PIN Food and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% Food shellings
% r=zeros(71,1);
% for ii=1:71
% T=S(ii:end);
% for jj=(ii-1):-1:1
% T=T+S(jj:end-(ii-jj));
% end
% r(ii)=corr(T,Fs(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shellings','Fontsize',22)
% box off;
% title('PIN Food and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
% r(ii)=corr(S(1:(end-(ii-1))),Fs(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of shellings','Fontsize',22)
% box off;
% title('PIN Food and shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% %% WaSH conflict gov. level
% r=zeros(71,1);
% for ii=1:71
% T=Cg(:,ii:end);
% for jj=(ii-1):-1:1
% T=T+Cg(:,jj:end-(ii-jj));
% end
% WW=Wsg(:,ii:end);
% r(ii)=corr(T(:),WW(:));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of conflict','Fontsize',22)
% box off;
% title('PIN WaSH and Conflict (Gov)');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
%     T=Cg(:,1:(end-(ii-1)));
%     
% WW=Wsg(:,ii:end);
% r(ii)=corr(T(:),WW(:));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of conflict','Fontsize',22)
% box off;
% title('PIN WaSH and Conflict (Gov)');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% 
% %% WaSH shellings gov. level
% r=zeros(71,1);
% for ii=1:71
% T=Sg(:,ii:end);
% for jj=(ii-1):-1:1
% T=T+Sg(:,jj:end-(ii-jj));
% end
% WW=Wsg(:,ii:end);
% r(ii)=corr(T(:),WW(:));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shelling','Fontsize',22)
% box off;
% title('PIN WaSH and shellings (Gov)');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% r=zeros(71,1);
% for ii=1:71
%     T=Sg(:,1:(end-(ii-1)));
%     
% WW=Wsg(:,ii:end);
% r(ii)=corr(T(:),WW(:));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Lag weeks of shelling','Fontsize',22)
% box off;
% title('PIN WaSH and shelling (Gov)');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% clear;
% % %% Conflict and Diesel
% % [WI,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,VT1,VT2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
% % 
% % D=Dieselt(GNZI,:);
% % W=Wheatt(GNZI,:);
% % C=Ctv(GNZI,:);
% % S=Mt(GNZI,:);
% % 
% % C=sum(C,1)';
% % S=sum(S,1)';
% % D=mean(D,1)';
% % W=mean(W,1)';
% % 
% % %% Diesel conflict
% % r=zeros(71,1);
% % for ii=1:71
% % T=C(ii:end);
% % for jj=(ii-1):-1:1
% % T=T+C(jj:end-(ii-jj));
% % end
% % r(ii)=corr(T,D(ii:end));
% % end
% % 
% % 
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% % ylabel('Correlation','Fontsize',22);
% % xlabel('Cumulative weeks of conflict','Fontsize',22)
% % box off;
% % title('Diesel and Conflict');
% % set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% % 
% % r=zeros(71,1);
% % for ii=1:71
% % r(ii)=corr(C(1:(end-(ii-1))),D(ii:end));
% % end
% % 
% % 
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% % ylabel('Correlation','Fontsize',22);
% % xlabel('Lag weeks of conflict','Fontsize',22)
% % box off;
% % title('Diesel and Conflict');
% % set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% % 
% % 
% % %% Diesel and shellings
% % 
% % r=zeros(71,1);
% % for ii=1:71
% % T=S(ii:end);
% % for jj=(ii-1):-1:1
% % T=T+S(jj:end-(ii-jj));
% % end
% % r(ii)=corr(T,D(ii:end));
% % end
% % 
% % 
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% % ylabel('Correlation','Fontsize',22);
% % xlabel('Cumulative weeks of shellings','Fontsize',22)
% % box off;
% % title('Diesel and Shellings');
% % set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% % 
% % r=zeros(71,1);
% % for ii=1:71
% % r(ii)=corr(S(1:(end-(ii-1))),D(ii:end));
% % end
% % 
% % 
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% % ylabel('Correlation','Fontsize',22);
% % xlabel('Lag weeks of shellings','Fontsize',22)
% % box off;
% % title('Diesel and shellings');
% % set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% % 
% % %% wheat and conflict
% % r=zeros(71,1);
% % for ii=1:71
% % T=C(ii:end);
% % for jj=(ii-1):-1:1
% % T=T+C(jj:end-(ii-jj));
% % end
% % r(ii)=corr(T,W(ii:end));
% % end
% % 
% % 
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% % bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% % ylabel('Correlation','Fontsize',22);
% % xlabel('Cumulative weeks of conflict','Fontsize',22)
% % box off;
% % title('Wheat and Conflict');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');
% 
% %% wheat and shellings
% 
% r=zeros(71,1);
% for ii=1:71
% T=S(ii:end);
% for jj=(ii-1):-1:1
% T=T+S(jj:end-(ii-jj));
% end
% r(ii)=corr(T,W(ii:end));
% end
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% bar([1:71],r,'Facecolor','k','Facealpha',0.35,'LineStyle','none')
% ylabel('Correlation','Fontsize',22);
% xlabel('Cumulative weeks of shellings','Fontsize',22)
% box off;
% title('Wheat and Shellings');
% set(gca,'linewidth',2,'tickdir','out','Fontsize',18,'Xminortick','on');