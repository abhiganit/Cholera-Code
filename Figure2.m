close all;
clear;
FC=[hex2rgb('#F5BE41');
    hex2rgb('#a6bddb');];
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_Gov_Incidence.mat')
IData=IData';
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.065,0.62,0.93,0.35]);
NW=length(WI(1,:));
bar([(1+maxtau):NW],sum(MI),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on
scatter([1:NW],sum(IData),40,'k','filled'); 
box off;
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTick',[0:5000:55000],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
hy=ylabel('Suspected cholera cases','Fontsize',18);
xlabel('Week reported','Fontsize',18);
IW=[1 21; 22 74; 75 121; 122 149];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',20);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',20);
ylim([0 55005]);
text(1,52005,'Model estimate','Fontsize',18,'Color',FC(1,:));
text(1,48411,'Reported','Fontsize',18,'Color','k');
text(hy.Extent(1),56020,'A','Fontsize',32,'FontWeight','bold');

print(gcf,['Figure2A.png'],'-dpng','-r600');

figure('units','normalized','outerposition',[0 0 1 1]);

subplot('Position',[0.065,0.62,0.46,0.35]);
dW=10;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
NW=length(WI(1,:));
bar([(1+maxtau):NW],MI(2,:),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on
scatter([1:NW],IData(2,:),20,'k','filled'); 
box off;
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTick',[0:1000:8000],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
hy=ylabel('Suspected cholera cases','Fontsize',18);
ylim([0 8000]);
text(1.5,7700,'\it{\bf{Aden}}','Fontsize',18);


text(1,7042,'Model estimate','Fontsize',18,'Color',FC(1,:));
text(1,6520,'Model validation','Fontsize',18,'Color',FC(2,:));
text(1,5998,'Reported','Fontsize',18,'Color','k');

            xlabel('Week reported','Fontsize',18);
 
text(hy.Extent(1),8000,'B','Fontsize',32,'FontWeight','bold');
subplot('Position',[0.535,0.62,0.46,0.35]);

NW=length(WI(1,:));
bar([(1+maxtau):NW],MI(9,:),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on
scatter([1:NW],IData(9,:),20,'k','filled'); 
box off;
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTickLabel','','YTick',[0:1000:8000],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
ylim([0 8000]);
text(1.5,7700,'\it{\bf{Amanat Al Asimah}}','Fontsize',18);

            xlabel('Week reported','Fontsize',18);

%% District panels


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



[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;
[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_District_Incidence.mat')
IData=IData';
load('PopulationSize_DistrictYemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019


endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end); % Incidence data for the districts starts at may 1 2017


MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
%% Aden
cc=[4 7];
for yy=1:1
    for xx=1:2
        subplot('Position',[0.095+0.225.*(xx-1),0.28-0.175*(yy-1),0.205,0.20]);
        NW=length(WI(1,:));
        bar([(1+maxtau):NW],MI(3+length(fS)+cc(xx),:),'Facecolor',FC(2,:),'LineStyle','none','Facealpha',1); hold on
        scatter([1:NW],IData(3+length(fS)+cc(xx),:),10,'k','filled'); 
        box off;
        dW=20;
        XTL=datestr([endDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
        xlim([0.5 length(WI(1,:))+0.5]);
%         if(yy~=4)
%             if(xx==1)
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:100:600],'YMinortick','on');
%                 
%             ylabel({'Suscpeted cholera cases'},'Fontsize',12);
%             else
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:100:600],'YMinortick','on','YTickLabel','');
%             end
%         else
            
            xlabel('Week reported','Fontsize',12);
            if(xx==1)
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:100:600],'YMinortick','on');
                
                ylabel({'Suscpeted cholera cases'},'Fontsize',12);
                text(-47.33,600,'C','Fontsize',32,'FontWeight','bold');
            else
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:100:600],'YMinortick','on','YTickLabel','');
            end
%         end
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        xtickangle(90);
        ylim([0 600]);
        text(1.5,580,['\it{\bf{' SD(fA(cc(xx))).ADM2_EN '}}'],'Fontsize',14);
        if(xx==1)
            text(79,548,'Model validation','Fontsize',14,'Color',FC(2,:));
            text(79,501,'Reported','Fontsize',14,'Color','k'); 
        end
    end
end

%% Sana'a City
cc=[2 7];
for yy=1:1
    for xx=1:2
        subplot('Position',[0.565+0.225.*(xx-1),0.28-0.175*(yy-1),0.205,0.20]);
        NW=length(WI(1,:));
        bar([(1+maxtau):NW],MI(3+cc(xx),:),'Facecolor',FC(2,:),'LineStyle','none','Facealpha',1); hold on
        scatter([1:NW],IData(3+cc(xx),:),10,'k','filled'); 
        box off;
        dW=20;
        XTL=datestr([endDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
        xlim([0.5 length(WI(1,:))+0.5]);
%         if(yy~=5)
%             if(xx==1)
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on');
%                 
%             ylabel({'Suscpeted cholera cases'},'Fontsize',12);
%             else
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on','YTickLabel','');
%             end
%         else
            xlabel('Week reported','Fontsize',12);
            if(xx==1)
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on');
                
             ylabel({'Suscpeted cholera cases'},'Fontsize',12);
            else
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on','YTickLabel','');
            end
%         end
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        xtickangle(90);
        ylim([0 1750]);
        text(1.5,1650,['\it{\bf{' SD(fS(cc(xx))).ADM2_EN '}}'],'Fontsize',14);
        
    end
end

print(gcf,['Figure2B.png'],'-dpng','-r600');
%% Hodeidah City
figure('units','normalized','outerposition',[0 0 1 1]);
dW=10;
subplot('Position',[0.535,0.62,0.46,0.35]);

NW=length(WI(1,:));
bar([(1+maxtau):NW],MI(3+length(fS)+length(fA)+1,:),'Facecolor',FC(2,:),'LineStyle','none','Facealpha',1); hold on
scatter([1:NW],IData(3+length(fS)+length(fA)+1,:),20,'k','filled'); 
box off;
xlim([0.5 length(WI(1,:))+0.5]);
set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on','YTickLabel','','YTick',[0:1000:8000],'YMinortick','on');
% Sets the y-axis to not have 10^n
ax=gca; % finds the current axis
ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
xtickangle(90);
ylim([0 8000]);
text(1.5,7700,['\it{\bf{Hodeidah City}}'],'Fontsize',18);
xlabel('Week reported','Fontsize',18);


%% Hodeidah City
cc=[1 2];
VGI=[29 31 71];
for yy=1:1
    for xx=1:2
        subplot('Position',[0.565+0.225.*(xx-1),0.28-0.175*(yy-1),0.205,0.20]);
        NW=length(WI(1,:));
        bar([(1+maxtau):NW],MI(cc(xx),:),'Facecolor',FC(2,:),'LineStyle','none','Facealpha',1); hold on
        scatter([1:NW],IData(cc(xx),:),10,'k','filled'); 
        box off;
        dW=20;
        XTL=datestr([endDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
        xlim([0.5 length(WI(1,:))+0.5]);
%         if(yy~=5)
%             if(xx==1)
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on');
%                 
%             ylabel({'Suscpeted cholera cases'},'Fontsize',12);
%             else
%                 set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',12,'Xminortick','on','YTick',[0:250:1750],'YMinortick','on','YTickLabel','');
%             end
%         else
            xlabel('Week reported','Fontsize',12);
            if(xx==1)
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:500:2500],'YMinortick','on');
                
             ylabel({'Suscpeted cholera cases'},'Fontsize',12);
            else
                set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',12,'Xminortick','on','YTick',[0:500:2500],'YMinortick','on','YTickLabel','');
            end
%         end
        % Sets the y-axis to not have 10^n
        ax=gca; % finds the current axis
        ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
        xtickangle(90);
        ylim([0 2500]);
        text(1.5,2400,['\it{\bf{' SD(VGI(cc(xx))).ADM2_EN '}}'],'Fontsize',14);
        
    end
end
print(gcf,['Figure2C.png'],'-dpng','-r600');