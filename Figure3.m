% Read table of past fitsclose all;
close all;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,VT1,VT2,GNZI,GV,maxtau] = LoadYemenData; % Load the data used to construct the figure
PDS=0.85; % Specify the size of the trainign set used
atest=0; % The selction criteria
load(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF,mln,a,KV,dV]=RetParameterPS(parv(end,:),XUv(end,:));

% Run the logistic model
[Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a,VT1(GNZI,1:length(WI(1,:))),VT2(GNZI,1:length(WI(1,:))),KV,dV);
% Compute vaccination
dV1=ImpactAttack(VT1(GNZI,1:length(WI(1,:))),0,dV,2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(VT2(GNZI,1:length(WI(1,:))),0,dV,2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV(1),dV2,KV(2));
% Compute conariates with and without WaSH and Food security
[XRemoveFS] = CalcCovariates(WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),WPIN(GNZI,1:length(WI(1,:))),0.*FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a);
[XRemoveWaSH] = CalcCovariates(WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI,1:length(WI(1,:))),RC(GNZI),H(GNZI,1:length(WI(1,:))),0.*WPIN(GNZI,1:length(WI(1,:))),FPIN(GNZI,1:length(WI(1,:))),Mt(GNZI,1:length(WI(1,:))),Wheatt(GNZI,1:length(WI(1,:))),Dieselt(GNZI,1:length(WI(1,:))),mln,a);

% load the populatino szie to determine the weekly incidence
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
% Compute the incidence 
MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
 
Inc=1:6;
FW=[1 4 5 6]; % Compentens that are dependent on the food and wash
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
     CIncFS=CIncFS+(1-EOVC).*beta(FW(ii)).*squeeze(XRemoveWaSH(FW(ii),:,:)).*PopS(GNZI,maxtau+1:end)./10000;
     CIncWaSH=CIncWaSH+(1-EOVC).*beta(FW(ii)).*squeeze(XRemoveFS(FW(ii),:,:)).*PopS(GNZI,maxtau+1:end)./10000;
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

IndW=[1 21; 22 74; 75 121; 122 149]; % Index of wave for the data used in the regression model
WW=zeros(4,length(GNZI),2);
for ww=1:4
    WW(ww,:,:)=[(sum(CIncFS(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2)) (sum(CIncWaSH(:,IndW(ww,1):IndW(ww,2)),2)./sum(MI(:,IndW(ww,1):IndW(ww,2)),2))];
end

%% Plot the data

ColorM=[0 0.6 0.6; % WaSH
        hex2rgb('#2E4600'); ];% Food Security 


IW=7.*(([1; 22 ; 75 ; 122; 150]-1)+maxtau); % The 150 is the start of the week we do not have data for and we are subtracting a week for the index of the week as the index zero is Oct 3, 2016
IW=[IW(1) IW(2)-1 IW(2) IW(3)-1 IW(3) IW(4)-1 IW(4) IW(5)-1];
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+IW],'mmm.dd,yyyy');

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

XGL={S(GNZI).ADM1_EN};



% Governorate
figure('units','normalized','outerposition',[0 0 1 1]);
% First wave
subplot('Position',[0.045,0.48,0.45,0.45]);

b=bar([1:length(GNZI)],squeeze(WW(1,:,:)),'LineStyle','none');
set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
xtickangle(45);
box off;
xlabel('Governorate','Fontsize',18);
title([XTL(1,:) ' to ' XTL(2,:)],'Fontsize',18);
for ii=1:length(ColorM(:,1))
    b(ii).FaceColor = 'flat';
    b(ii).CData = ColorM(ii,:);
end
ylim([0 1]);
yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);
ax=gca;
tickl=ax.TickLength;
text(yh.Extent(1),max(ylim)*1.05,'A','Fontsize',32,'FontWeight','bold');
% Second wave
subplot('Position',[0.545,0.48,0.45,0.45]);

b=bar([1:length(GNZI)],squeeze(WW(2,:,:)),'LineStyle','none');
set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
xtickangle(45);
box off;
xlabel('Governorate','Fontsize',18);
title([XTL(3,:) ' to ' XTL(4,:)],'Fontsize',18);
for ii=1:length(ColorM(:,1))
    b(ii).FaceColor = 'flat';
    b(ii).CData = ColorM(ii,:);
end

ylim([0 1]);
yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);

text(yh.Extent(1),max(ylim)*1.05,'B','Fontsize',32,'FontWeight','bold');

figure('units','normalized','outerposition',[0 0 1 1]);
% Third wave
subplot('Position',[0.045,0.48,0.45,0.45]);

b=bar([1:length(GNZI)],squeeze(WW(3,:,:)),'LineStyle','none');
set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
xtickangle(45);
box off;
xlabel('Governorate','Fontsize',18);
title([XTL(5,:) ' to ' XTL(6,:)],'Fontsize',18);
for ii=1:length(ColorM(:,1))
    b(ii).FaceColor = 'flat';
    b(ii).CData = ColorM(ii,:);
end

ylim([0 1]);
yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);

text(yh.Extent(1),max(ylim)*1.05,'C','Fontsize',32,'FontWeight','bold');
% Fourth wave
subplot('Position',[0.545,0.48,0.45,0.45]);

b=bar([1:length(GNZI)],squeeze(WW(4,:,:)),'LineStyle','none');
set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',XGL,'Fontsize',16,'Yminortick','on');
xtickangle(45);
box off;
xlabel('Governorate','Fontsize',18);
title([XTL(7,:) ' to ' XTL(8,:)],'Fontsize',18);
for ii=1:length(ColorM(:,1))
    b(ii).FaceColor = 'flat';
    b(ii).CData = ColorM(ii,:);
end

ylim([0 1]);
yh=ylabel('Contribution to suscpected cholera cases','Fontsize',18);

    text(yh.Extent(1),max(ylim)*1.05,'D','Fontsize',32,'FontWeight','bold');
% % National
% figure('units','normalized','outerposition',[0 0 1 1]);
% 
% subplot('Position',[0.045,0.48,0.955,0.45]);
% 
% 
% b=bar([1:4],rand(4,3),'LineStyle','none');
% ylim([0 1]);
% set(gca,'linewidth',2,'tickdir','out','XTick',[1:length(GNZI)],'XTickLabel',{'First','Second','Third','Fourth'},'Fontsize',16,'Yminortick','on','TickLength',tickl.*0.45/0.955);
% 
% box off;
% xlabel('Epidemic wave','Fontsize',18);
% ylabel('Contribution to suscpected cholera cases','Fontsize',18);
% for ii=1:length(ColorM(:,1))
%     b(ii).FaceColor = 'flat';
%     b(ii).CData = ColorM(ii,:);
% end