close all;
clear;
FC=[hex2rgb('#F5BE41');
    hex2rgb('#a6bddb');];
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDataVal;
NW=166; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
NW1=153;
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);


[Yt,~]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
load('Yemen_Gov_Incidence_Val.mat')
IData=IData';
IData=IData(GNZI,:);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=166-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.065,0.62,0.93,0.35]);
NW=length(WI(1,:));
b=bar([(1+maxtau):NW],sum(MI),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on

b.FaceColor = 'flat';
for mm=((NW1-maxtau)+1:162)
    b.CData(mm,:) = FC(2,:);
end

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
IW=[1 21; 22 74; 75 121; 122 162];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',18);
ylim([0 55005]);
print(gcf,['TemporalVal_Yemen.png'],'-dpng','-r600');

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

Location={S(GNZI).ADM1_EN}';
llb=[0 500 2000 5000];
uub=[500 2000 5000 8000];
dx=[50 250 500 1000];
CLO=struct('N',{'Low','ModLow','ModHigh','High'});
for ll=1:4
    findx=find(max(MI,[],2)<uub(ll));
    gindx=find(max(MI(findx,:),[],2)>=llb(ll));
    findx=findx(gindx);
    gmax=ceil(length(findx)/4);
    for gg=1:gmax
        figure('units','normalized','outerposition',[0 0 1 1]);
        for ss=1:4        
            if(ss+4.*(gg-1)<=length(findx))
                subplot('Position',[0.065*(rem(ss,2))+0.535*(1-rem(ss,2)),0.6*(1-floor(ss/3))+0.19*floor(ss/3),0.46,0.35]);
                GGG=findx(ss+4.*(gg-1));    
                dW=10;
                XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
                NW=length(WI(1,:));
                b=bar([(1+maxtau):NW],MI(GGG,:),'Facecolor',FC(1,:),'LineStyle','none','Facealpha',1); hold on

                b.FaceColor = 'flat';
                for mm=((NW1-maxtau)+1:162)
                    b.CData(mm,:) = FC(2,:);
                end
                scatter([1:NW],IData(GGG,:),20,'k','filled'); 
                box off;
                xlim([0.5 length(WI(1,:))+0.5]);
                if((gg==gmax)&&(ss>2))||((ss+4.*(gg-1))==length(findx))||((ss+4.*(gg-1))==(length(findx)-1))||((gg==(gmax-1))&&((ss+4.*(gg-1))==(length(findx)-2)))
                    if(rem(ss,2)==1)
                        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','off','YTick',[0:dx(ll):uub(ll)],'YMinortick','on');
                        xlabel('Week reported','Fontsize',18);
                        hy=ylabel('Suspected cholera cases','Fontsize',18);
                    else                
                        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','off','YTick',[0:dx(ll):uub(ll)],'YMinortick','on','YTickLabel','');
                        xlabel('Week reported','Fontsize',18);
                        hy=ylabel('','Fontsize',18);
                    end
                else
                    if(rem(ss,2)==1)
                        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',16,'Xminortick','off','YTick',[0:dx(ll):uub(ll)],'YMinortick','on');
                        hy=ylabel('Suspected cholera cases','Fontsize',18);
                    else                
                        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel','','Fontsize',16,'Xminortick','off','YTick',[0:dx(ll):uub(ll)],'YMinortick','on','YTickLabel','');
                        hy=ylabel('','Fontsize',18);
                    end
                end
                % Sets the y-axis to not have 10^n
                ax=gca; % finds the current axis
                ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
                xtickangle(90);
                ylim([0 uub(ll)]);
                text(1.5,7700/8000*uub(ll),Location{GGG},'Fontsize',18);
                text(hy.Extent(1),8750/8000*uub(ll),char(64+ss+4.*(gg-1)),'Fontsize',32,'FontWeight','bold');
            end
        end
        print(gcf,['TemporalVal_' CLO(ll).N '_' num2str(gg) '.png'],'-dpng','-r600');
    end
end
