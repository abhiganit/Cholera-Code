close all;
clc;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
RC=RC(GNZI);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);
beta(13:16)
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
XC2=zeros(size(squeeze(X(1,:,:))));
for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
    for gg=1:21
        XC2(gg,:)=pchip([1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))],squeeze(XC(ii-maxtau*(mmt-1),gg,:)),[1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))]-1);
    end
    tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*(bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
    tempmat2=tempmat2+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
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
WWT=[WW(:,1:6) sum(WW(:,2:3),2)];
[WS, indexs]=sortrows(WWT,7);
WS=WS(:,1:6);
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
h=ylabel('Governorate','Fontsize',18);
text(h.Extent(1),max(ylim),'A','Fontsize',32,'Fontweight','bold');
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

%% District
clear;
clc;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
GNZI=[4:13];
NW=123; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
RC=RC(GNZI);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;
[Yt,X]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);

load('PopulationSize_DistrictYemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019


endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end); % Incidence data for the districts starts at may 1 2017

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
load('DieselrepresentedthroughConflictShellings.mat','bd');
load('DieselrepresentedthroughConflictShellings_District.mat','XC','XS');
XC=XC(:,GNZI,:);
XS=XS(:,GNZI,:);
mmt=4;
tempmat=zeros(size(squeeze(X(1,:,:))));
tempmat2=zeros(size(squeeze(X(1,:,:))));
XC2=zeros(size(squeeze(X(1,:,:))));
for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
    for gg=1:length(GNZI)
        XC2(gg,:)=pchip([1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))],squeeze(XC(ii-maxtau*(mmt-1),gg,:)),[1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))]-1);
    end
    tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*(bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
    tempmat2=tempmat2+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
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
    for jj=1:length(GNZI)
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


SD=SD([29 31 71 fS' fA']);

XGL={SD(GNZI).ADM2_EN};
WWT=[WW(:,1:6) sum(WW(:,2:3),2)];
[WS, indexs]=sortrows(WWT,7);
WS=WS(:,1:6);
subplot('Position',[0.6,0.59,0.347689075630252,0.32]);
b=barh([1:length(GNZI)],WS,'stacked','LineStyle','none'); 
for ii=1:6
   b(ii).FaceColor=ColorM(ii,:); 
end
XGL2=XGL;
for ii=1:length(GNZI)
    XGL2{ii}=XGL{indexs(ii)};
end
hold on
dx1=linspace(0,0.45,101);
dx2=dx1+(dx1(2)-dx1(1))./2;
dx2=dx2(1:end-1);
for mm=1:length(GNZI)
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
h=ylabel('Districts in Amanat Al Asimah','Fontsize',18);
text(h.Extent(1),max(ylim),'B','Fontsize',32,'Fontweight','bold');
xlim([0 0.45]);
ylim([0.5 length(GNZI)+.5]);
box off;
set(gca,'LineWidth',2,'tickdir','out','YTick',[1:length(GNZI)],'YTickLabel',XGL2,'Fontsize',16);
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
line((CI(GNZI(indexs),end)),[1:length(GNZI)],'Parent',ax2,'Color',[0.5 0.5 0.5],'LineWidth',2)
ylim([0.5 length(GNZI)+.5]);
xlim([10^2 10^5])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:length(GNZI)],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence','Fontsize',18,'Color',[0.5 0.5 0.5]);

%% District (Aden)
clear;
clc;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat')

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
GNZI=[14:21];
NW=123; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);
RC=RC(GNZI);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7)-1;

n=n-TruncV;
n(n<1)=1;
[Yt,X]= LogisticModel(beta,tA(GNZI,:),DB,DA,Ctv(GNZI,:),K,n,tau,maxtau,CF,WPIN(GNZI,:),FPIN(GNZI,:),Mt(GNZI,:),Wheatt(GNZI,:),Dieselt(GNZI,:),KP,V1(GNZI,:),V2(GNZI,:),KV,dV,Rtv(GNZI,:),RF,r0,WI(GNZI,:),PopS(GNZI,:),CI(GNZI,:),DAR,w);
dV1=ImpactAttack(V1(GNZI,:)-V2(GNZI,:),0,dV(1),2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(V2(GNZI,:),0,dV(2),2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV,dV2,KV);

load('PopulationSize_DistrictYemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019


endDateofSim = datenum('5-01-2017');% End date
TruncV=ceil((1+endDateofSim-startDateofSim)./7);

% External effect due to IDP
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; % population size to feed into the IDPt calculation
PopS=PopS(:,TruncV:end); % Incidence data for the districts starts at may 1 2017

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
load('DieselrepresentedthroughConflictShellings.mat','bd');
load('DieselrepresentedthroughConflictShellings_District.mat','XC','XS');
XC=XC(:,GNZI,:);
XS=XS(:,GNZI,:);
mmt=4;
tempmat=zeros(size(squeeze(X(1,:,:))));
tempmat2=zeros(size(squeeze(X(1,:,:))));
XC2=zeros(size(squeeze(X(1,:,:))));
for ii=(maxtau*(mmt-1)+1):(mmt.*maxtau)
    for gg=1:length(GNZI)
        XC2(gg,:)=pchip([1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))],squeeze(XC(ii-maxtau*(mmt-1),gg,:)),[1:length(squeeze(XC(ii-maxtau*(mmt-1),gg,:)))]-1);
    end
    tempmat=tempmat+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*(bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
    tempmat2=tempmat2+(1-EOVC).*(beta(ii).*squeeze(X(ii,:,:))).*PopS(GNZI,maxtau+1:end)./10000.*bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))./(bd(1)+bd(2).*squeeze(XC(ii-maxtau*(mmt-1),:,:))+bd(3).*squeeze(XS(ii-maxtau*(mmt-1),:,:))+bd(4).*squeeze(XC2));
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
    for jj=1:length(GNZI)
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


SD=SD([29 31 71 fS' fA']);

XGL={SD(GNZI).ADM2_EN};
WWT=[WW(:,1:6) sum(WW(:,2:3),2)];
[WS, indexs]=sortrows(WWT,7);
WS=WS(:,1:6);
subplot('Position',[0.6,0.109422492401216,0.347689075630252,0.32]);
b=barh([1:length(GNZI)],WS,'stacked','LineStyle','none'); 
for ii=1:6
   b(ii).FaceColor=ColorM(ii,:); 
end
XGL2=XGL;
for ii=1:length(GNZI)
    XGL2{ii}=XGL{indexs(ii)};
end
hold on
dx1=linspace(0,0.45,101);
dx2=dx1+(dx1(2)-dx1(1))./2;
dx2=dx2(1:end-1);
for mm=1:length(GNZI)
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
h=ylabel('Districts in Aden','Fontsize',18);
text(h.Extent(1),max(ylim),'C','Fontsize',32,'Fontweight','bold');
xlim([0 0.45]);
ylim([0.5 length(GNZI)+.5]);
box off;
set(gca,'LineWidth',2,'tickdir','out','YTick',[1:length(GNZI)],'YTickLabel',XGL2,'Fontsize',16);
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
line((CI(GNZI(indexs),end)),[1:length(GNZI)],'Parent',ax2,'Color',[0.5 0.5 0.5],'LineWidth',2)
ylim([0.5 length(GNZI)+.5]);
xlim([10^2 10^5])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:length(GNZI)],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence','Fontsize',18,'Color',[0.5 0.5 0.5]);
