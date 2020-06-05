close all;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
RC=RC(GNZI);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);

[Yt,X]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);
load('PopulationSize_Yemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation

MI=(Yt./(10000)).*PopS(GNZI,maxtau+1:end);
 CCR=cell(6,1);
 for mm=1:6
     tempmat=zeros(size(squeeze(X(1,:,:))));
    for ii=(maxtau*(mm-1)+1):(mm.*maxtau)
        tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000;
    end
    CCR{mm}=tempmat;
 end
%% Conflict indirect effect
load('DieselrepresentedthroughConflictShellings.mat','bd','XC','XS');
mmt=4;
tempmat=zeros(size(squeeze(X(1,:,:))));
tempmat2=zeros(size(squeeze(X(1,:,:))));
for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
    tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:)));
    tempmat2=tempmat2+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:)));
end
CCR{2}=CCR{2}+tempmat;
CCR{3}=CCR{3}+tempmat2;
CCR{4}=CCR{4}-tempmat-tempmat2;

IndW=[1 21; 22 74; 75 121; 122 149]; % Index of wave for the data used in the regression model
WW=zeros(length(GNZI),7);
WW2=zeros(length(GNZI),6);
for mm=1:6
    temp=((CCR{mm}))./(MI);
    temp3=((CCR{mm}))./(CCR{1}+CCR{2}+CCR{3}+CCR{4}+CCR{5}+CCR{6});
    for jj=1:21
        temp2=temp(jj,MI(jj,:)>0);
        WW(jj,mm) = mean(temp2);
        temp4=temp3(jj,MI(jj,:)>0);
        WW2(jj,mm) = mean(temp4);
    end
end  

mm=7;
    WW(:,mm)=1-sum(WW(:,1:6),2);



explode = [0,1,1,0,0,0];
ColorM=[[221,28,119]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        ]; 
    labels = {'Targeted','Conflict','Shellings','Diesel','Wheat','Rainfall'};
    
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

XGL={S(GNZI).ADM1_EN};
[WS, indexs]=sortrows(WW(:,1:6),2);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.117121848739496,0.109422492401216,0.347689075630252,0.8]);
b=barh([1:21],WS,'stacked','LineStyle','none'); 
for ii=1:6
   b(ii).FaceColor=ColorM(ii,:); 
end
XGL2=XGL;
for ii=1:21
    XGL2{ii}=XGL{indexs(ii)};
end
hold on
dx1=linspace(0,0.45,101);
dx2=dx1+(dx1(2)-dx1(1))./2;
dx2=dx2(1:end-1);
for mm=1:21
    if(RC(indexs(mm))==1)
        for ii=1:5
            if((ii==1)||(ii==3)||(ii==5))
               scatter(dx1(dx1<sum(WS(mm,1:6))),(0.15.*(ii-3)+ mm).*ones(size(dx1(dx1<sum(WS(mm,1:6))))),5,'k','filled');
            else
                scatter(dx2(dx2<sum(WS(mm,1:6))),(0.15.*(ii-3)+ mm).*ones(size(dx2(dx2<sum(WS(mm,1:6))))),5,'k','filled');
            end
        end
    end
end
xlabel('Average relative contribution per week','Fontsize',18);
ylabel('Governorate','Fontsize',18);
xlim([0 0.45]);
ylim([0.5 21.5]);
box off;
set(gca,'LineWidth',2,'tickdir','out','YTick',[1:21],'YTickLabel',XGL2,'Fontsize',16);
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
line((CI(GNZI(indexs),end)),[1:21],'Parent',ax2,'Color',[0.5 0.5 0.5],'LineWidth',2)
ylim([0.5 21.5]);
xlim([10^2 10^6])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:21],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence','Fontsize',18,'Color',[0.5 0.5 0.5]);

idx = kmeans([WS CI(GNZI(indexs))],2);
test=0;
% for kk=1:4
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     if(kk<4)
%         for nn=1:2
%             for jj=1:3    
%                 subplot('Position',[0.14+(jj-1)*0.29,0.5-0.45.*(nn-1),0.27,0.4]);
%                 p=pie3(WW2(jj+3.*(nn-1)+6.*(kk-1),:),explode,labels);
%                 set(p(4:4:end),'FontSize',14);
%                 colormap(ColorM)
%                 TT=[XGL(jj+3.*(nn-1)+6.*(kk-1))];
%                 title({TT{1},[ '(' num2str(round(100.*(1-WW(jj+3.*(nn-1)+6.*(kk-1),7)),1)) '%)']},'Fontsize',18)
%             end
%         end
%     else
%         for jj=1:3    
%             subplot('Position',[0.14+(jj-1)*0.29,0.5-0.45.*(nn-1),0.27,0.4]);
%             p=pie3(WW2(jj+6.*(kk-1),:),explode,labels);
%             set(p(4:4:end),'FontSize',14);
%             colormap(ColorM)
%             TT=[XGL(jj+6.*(kk-1))];
%             title({TT{1},[ '(' num2str(round(100.*(1-WW(jj+6.*(kk-1),7)),1)) '%)']},'Fontsize',18)
%             if(jj==2)
%                 legend({'Conflict','Shellings','Diesel','Rainfall'},'Location','southoutside','Orientation','horizontal','Fontsize',18)
%                 legend boxoff;
%             end
%         end
%     end
%     print(gcf,['FigureRel' num2str(kk) '.png'],'-dpng','-r600');
% end
