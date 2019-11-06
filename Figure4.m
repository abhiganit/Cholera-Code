% Read table of past fitsclose all;
close all;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,GNZI,maxtau] = LoadYemenDistrictData; % Load the data used to construct the figure
PDS=0.65;
atest=0;
load(['ForwardSelectionNoRainNoConflict-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF,mln,a]=RetParameterPS(parv(end,:),XUv(end,:));

[Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a);

[XRemoveFS] = CalcCovariates(WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),0.*FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a);
[XRemoveWaSH] = CalcCovariates(WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),0.*WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a);

load('PopulationSize_DistrictYemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end);
MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
 
Inc=1:6;
FW=[1 4 5 6];
% Conf=7:10;
% Rain=11:12;
% CR=13:15;
% Com=16:18;
% CInc=zeros(size(squeeze(X(1,:,:))));
% CConflict=zeros(size(squeeze(X(1,:,:))));
% CRain=zeros(size(squeeze(X(1,:,:))));
% CCR=zeros(size(squeeze(X(1,:,:))));
% CCom=zeros(size(squeeze(X(1,:,:))));

CIncFS=zeros(size(squeeze(X(1,:,:))));
CIncWaSH=zeros(size(squeeze(X(1,:,:))));
for ii=1:length(FW)
     CIncFS=CIncFS+beta(FW(ii)).*squeeze(XRemoveWaSH(FW(ii),:,:)).*PopS(GNZI,maxtau+1:end)./10000;
     CIncWaSH=CIncWaSH+beta(FW(ii)).*squeeze(XRemoveFS(FW(ii),:,:)).*PopS(GNZI,maxtau+1:end)./10000;
end
% 
% for ii=1:length(Conf)
%     CConflict=CConflict+beta(Conf(ii)).*squeeze(X(Conf(ii),:,:));
% end
% 
% for ii=1:length(Rain)
%     CRain=CRain+beta(Rain(ii)).*squeeze(X(Rain(ii),:,:));
% end
% 
% for ii=1:length(CR)
%     CCR=CCR+beta(CR(ii)).*squeeze(X(CR(ii),:,:));
% end
% 
% for ii=1:length(Com)
%     CCom=CCom+beta(Com(ii)).*squeeze(X(Com(ii),:,:));
% end
%    
% PC=CCR./Yt;
% PC(Yt==0)=0;



IndW=[31+max(tau) 74; 75 121; 122 149]-(TruncV-1); % Index of wave for the data used in the regression model
WW=zeros(3,length(GNZI),2);
for ww=1:3
    WW(ww,:,:)=[(sum(CIncFS(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2)) (sum(CIncWaSH(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2))];
end

%% Plot the data

ColorM=[0 0.6 0.6; % WaSH
        hex2rgb('#2E4600'); ];% Food Security
SD = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen

fS=zeros(length(SD),1);
for ii=1:length(fS)
  fS(ii)=strcmp({'Amanat Al Asimah'},SD(ii).ADM1_EN); 
end
fS=find(fS==1);

fA=zeros(length(SD),1);
for ii=1:length(fA)
  fA(ii)=strcmp({'Aden'},SD(ii).ADM1_EN); 
end

fA=find(fA==1);



HodeidahCity=[29 31 71];

IW=7.*(([TruncV ; 75 ; 122; 150]-1)+maxtau); % The 150 is the start of the week we do not have data for and we are subtracting a week for the index of the week as the index zero is Oct 3, 2016 (31 is the 
IW=[IW(1) IW(2)-1 IW(2) IW(3)-1 IW(3) IW(4)-1 IW(4)];
dW=5;
XTL=datestr([startDateofSim+IW],'mmm.dd,yyyy');

%% Dates and indexs for the map plots
temp=[ 75 ; 122;]+maxtau; % Starting index of third and fourth wave relative to Oct. 3, 2016
 IdXMap=[1 temp(1)-(TruncV-1)-1; % need to take the index before the start of the next wave
       temp(1)-(TruncV-1) temp(2)-(TruncV-1)-1; % need to take the index before the start of the next wave
       temp(2)-(TruncV-1) length(WPIN(1,:))];
 SW=datestr([endDateofSim+7.*(IdXMap(:,1)-1)],'mmm.dd,yyyy');
 EW=datestr([endDateofSim+7.*([IdXMap(2:3,1);length(WPIN(1,:))+1]-1)-1],'mmm.dd,yyyy');
 TW=datestr([endDateofSim+7.*(IdXMap(:,2)-1)],'mmm.dd,yyyy');

% Aden
XGL={SD(fA).ADM2_EN};
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.04,0.59,0.955,0.4]);
rtemp=[];
for ii=1:3
   rtemp=[rtemp;squeeze(WW(ii,(3+length(fS)+1):(3+length(fS)+length(fA)),:))];
end
rr=cumsum(rtemp,2);
wp=0.2;
dwave=0.05;
dd=1;
for ww=1:3
       for ii=1:length(fA)
          for jj=1:length(ColorM(:,1))
              if(jj==1)
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[0 rr(ii+length(fA).*(ww-1),jj) rr(ii+length(fA).*(ww-1),jj) 0],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              else                
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[rr(ii+length(fA).*(ww-1),jj-1) rr(ii+length(fA).*(ww-1),jj) rr(ii+length(fA).*(ww-1),jj) rr(ii+length(fA).*(ww-1),jj-1)],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              end
          end
       end
end

set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(fA)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on','TickLength',[ 0.0047    0.0118]);
xtickangle(45);
box off;
xlabel('Districts in Aden','Fontsize',18);
ylabel('Suscpected cholera cases','Fontsize',18);
ylim([0 1]);
 xlim([0.6 length(fA)+0.4]);
 
   Malpha=zeros(3,length(fA));
 for ww=1:3
     for ii=1:length(fA)
         Malpha(ww,ii)=mean(a.*WPIN(3+length(fS)+ii,IdXMap(ww,1):IdXMap(ww,2))+(1-a).*FPIN(3+length(fS)+ii,IdXMap(ww,1):IdXMap(ww,2)));
     end
 end
 Malpha=(Malpha-min(Malpha(:)))./(max(Malpha(:))-min(Malpha(:)));
 for ww=1:3
     pp=subplot('Position',[0.04+(ww-1).*0.31,0.04,0.30,0.37]);
     for ii=1:length(fA)
            mapshow(SD(fA(ii)),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',Malpha(ww,ii)); hold on
     end
     xlim([44.1 45.11]);
     axis off;
     title([SW(ww,:) ' to ' EW(ww,:) ],'Fontsize',16);
     axes('Position',[-0.0225+(ww-1).*0.31,0.04 .2 .2])
     ii=3;
     mapshow(SD(fA(ii)),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',Malpha(ww,ii)); hold on
     box on;
     set(gca,'linewidth',2,'ticklength',[0 0],'YTickLabel','','XTickLabel','');
     xlim([43.36 43.45]);
     ylim([12.6 12.7]);
 end
 
 % Sana City
XGL={SD(fS).ADM2_EN};
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.04,0.59,0.955,0.4]);
rtemp=[];
for ii=1:3
   rtemp=[rtemp;squeeze(WW(ii,(3+1):(3+length(fS)),:))];
end
rr=cumsum(rtemp,2);
wp=0.2;
dwave=0.05;
dd=1;
for ww=1:3
       for ii=1:length(fS)
          for jj=1:length(ColorM(:,1))
              if(jj==1)
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[0 rr(ii+length(fS).*(ww-1),jj) rr(ii+length(fS).*(ww-1),jj) 0],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              else                
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[rr(ii+length(fS).*(ww-1),jj-1) rr(ii+length(fS).*(ww-1),jj) rr(ii+length(fS).*(ww-1),jj) rr(ii+length(fS).*(ww-1),jj-1)],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              end
          end
       end
end


set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(fS)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on','TickLength',[ 0.0047    0.0118]);
xtickangle(45);
box off;
xlabel('Districts in Sanaa city','Fontsize',18);
ylabel('Suscpected cholera cases','Fontsize',18);
ylim([0 1]);
 xlim([0.6 length(fS)+0.4]);

 
  Malpha=zeros(3,length(fS));
 for ww=1:3
     for ii=1:length(fS)
         Malpha(ww,ii)=mean(a.*WPIN(3+ii,IdXMap(ww,1):IdXMap(ww,2))+(1-a).*FPIN(3+ii,IdXMap(ww,1):IdXMap(ww,2)));
     end
 end
 Malpha=(Malpha-min(Malpha(:)))./(max(Malpha(:))-min(Malpha(:)));
 for ww=1:3
     subplot('Position',[0.04+(ww-1).*0.31,0.03,0.30,0.37]);
     for ii=1:length(fS)
        mapshow(SD(fS(ii)),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',Malpha(ww,ii)); hold on
     end
     axis off;
     title([SW(ww,:) ' to ' EW(ww,:) ],'Fontsize',16);
 end

 % Hodeidah City
XGL={SD(HodeidahCity).ADM2_EN};
figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.04,0.59,0.955,0.4]);
rtemp=[];
for ii=1:3
   rtemp=[rtemp;squeeze(WW(ii,1:3,:))];
end
rr=cumsum(rtemp,2);
wp=0.2;
dwave=0.05;
dd=1;
for ww=1:3
       for ii=1:length(HodeidahCity)
          for jj=1:length(ColorM(:,1))
              if(jj==1)
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[0 rr(ii+length(HodeidahCity).*(ww-1),jj) rr(ii+length(HodeidahCity).*(ww-1),jj) 0],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              else                
                p=patch(dd*ii+(wp+dwave)*(ww-2)+[-wp/2 -wp/2 wp/2 wp/2],[rr(ii+length(HodeidahCity).*(ww-1),jj-1) rr(ii+length(HodeidahCity).*(ww-1),jj) rr(ii+length(HodeidahCity).*(ww-1),jj) rr(ii+length(HodeidahCity).*(ww-1),jj-1)],ColorM(jj,:),'Facealpha',ww/3,'Edgealpha',0);
              end
          end
       end
end


set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(HodeidahCity)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on','TickLength',[ 0.0047    0.0118]);
xtickangle(45);
box off;
ylim([0 1]);
xlabel('Districts in Hodeidah City','Fontsize',18);
ylabel('Suscpected cholera cases','Fontsize',18);
 xlim([0.6 length(HodeidahCity)+0.4]);
 
 
 Malpha=zeros(3,length(HodeidahCity));
 for ww=1:3
     for ii=1:length(HodeidahCity)
         Malpha(ww,ii)=mean(a.*WPIN(ii,IdXMap(ww,1):IdXMap(ww,2))+(1-a).*FPIN(ii,IdXMap(ww,1):IdXMap(ww,2)));
     end
 end
 Malpha=(Malpha-min(Malpha(:)))./(max(Malpha(:))-min(Malpha(:)));
 for ww=1:3
     subplot('Position',[0.04+(ww-1).*0.31,0.04,0.30,0.37]);
     for ii=1:length(HodeidahCity)
        mapshow(SD(HodeidahCity(ii)),'FaceColor',hex2rgb('#2E4600'),'Edgecolor',[0 0 0],'LineWidth',2,'FaceAlpha',Malpha(ww,ii)); hold on
     end
     axis off;
     title([SW(ww,:) ' to ' EW(ww,:) ],'Fontsize',16);
 end
 