%% Plot Fuelt time line
close all;
clc;
clear;

load('Fuel_Yemen_Time.mat');

% Starts May 2014 to September 2019
DS=zeros(65,1);
cc=1;
for yy=2014:2019
    if(yy==2014)
        for mm=5:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    elseif(yy<2019)        
        for mm=1:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    else
       for mm=1:9
           DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
       end
    end
end

XTL=datestr(DS,'mm/dd/yy');
XTL([2:2:64],:)= ' ';

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

plot([1:65],FuelTime(:,1),'k','LineWidth',2); hold on;
plot([1:65],FuelTime(:,2),'color',hex2rgb('#CB0000'),'LineWidth',2);
plot([1:65],FuelTime(:,3),'color',hex2rgb('#F5BE41'),'LineWidth',2);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16,'XTickLabel',XTL,'XTick',[1:1:65],'YTick',[0:100:1500]);
xtickangle(45);
box off;
xlim([1 65]);
xlabel('Date','Fontsize',18);
ylabel('Price of fuel','Fontsize',18);
legend('Sana`a City','Hodeidah City','Aden');
legend box off;
%% Diesel
load('Diesel_Yemen_Time.mat');

% Starts May 2014 to September 2019
DS=zeros(65,1);
cc=1;
for yy=2014:2019
    if(yy==2014)
        for mm=5:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    elseif(yy<2019)        
        for mm=1:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    else
       for mm=1:9
           DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
       end
    end
end

XTL=datestr(DS,'mm/dd/yy');
XTL([2:2:64],:)= ' ';

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

plot([1:65],DieselTime(:,1),'k','LineWidth',2); hold on;
plot([1:65],DieselTime(:,2),'color',hex2rgb('#CB0000'),'LineWidth',2);
plot([1:65],DieselTime(:,3),'color',hex2rgb('#F5BE41'),'LineWidth',2);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16,'XTickLabel',XTL,'XTick',[1:1:65],'YTick',[0:100:1500]);
xtickangle(45);
box off;
xlim([1 65]);
xlabel('Date','Fontsize',18);
ylabel('Price of diesel','Fontsize',18);
legend('Sana`a City','Hodeidah City','Aden');
legend box off;

%% Diesel
load('Wheat_Yemen_Time.mat');

% Starts June 2016 to September 2019
DS=zeros(40,1);
cc=1;
for yy=2016:2019
    if(yy==2016)
        for mm=6:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    elseif(yy<2019)        
        for mm=1:12
            DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
        end
    else
       for mm=1:9
           DS(cc) = datenum(yy,mm,15);
            cc=cc+1;
       end
    end
end

XTL=datestr(DS,'mm/dd/yy');
XTL([2:2:40],:)= ' ';

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

plot([1:40],W(:,1),'k','LineWidth',2); hold on;
plot([1:40],W(:,2),'color',hex2rgb('#CB0000'),'LineWidth',2);
plot([1:40],W(:,3),'color',hex2rgb('#F5BE41'),'LineWidth',2);
set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16,'XTickLabel',XTL,'XTick',[1:1:40],'YTick',[0:25:300]);
xtickangle(45);
box off;
xlim([1 40]);
xlabel('Date','Fontsize',18);
ylabel('Price of wheat','Fontsize',18);
legend({'Sana`a City','Hodeidah City','Aden'},'Location','northwest');
legend box off;

%% Prior shellings

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt


load('PriorShelling_Yemen.mat');

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

ShellP(:,1)=ShellP(:,1)-min(ShellP(:,1))+1;

M=GLevelConflict(ShellP,SD,92);

load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,SD,153);

M=[M Mt];

plot([-91:153],[sum(M(3+[1:length(fS)],:),1)],'k','LineWidth',2); hold on;
plot([-91:153],[sum(M([1 2 3],:),1)],'color',hex2rgb('#CB0000'),'LineWidth',2)
plot([-91:153],[sum(M(3+length(fS)+[1:length(fA)],:),1)],'color',hex2rgb('#F5BE41'),'LineWidth',2)


startDateofSim = datenum('10-03-2016');% Start date
XTL=datestr([startDateofSim+7.*([-91:5:153])],'mm/dd/yy');

set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16,'XTickLabel',XTL,'XTick',([-91:5:153]));
xtickangle(45);
box off;
xlim([-91 153]);
xlabel('Date','Fontsize',18);
ylabel('Number of attacks','Fontsize',18);
legend('Sana`a City','Hodeidah City','Aden');
legend box off;

% Attack rate


figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,IDPt,GNZI,maxtau] = LoadYemenData;

[WId,Ctv,tA,Rtv,Mt,P,RC,H,WPINm,IDPt,GNZI,maxtau] = LoadYemenDistrictData;

startDateofSim = datenum('10-03-2016');% Start date
XTL=datestr([startDateofSim+7.*([0:5:152])],'mm/dd/yy');

plot([1:153],WI(9,:),'k','LineWidth',2); hold on;
plot([31:153],WId(22,:),'color',hex2rgb('#CB0000'),'LineWidth',2)
plot([1:153],WI(2,:),'color',hex2rgb('#F5BE41'),'LineWidth',2)

xtickangle(45);
box off;
xlim([1 153]);
ylim([0 110]);
xlabel('Date','Fontsize',18);
ylabel('Attack rate per 10,000','Fontsize',18);
legend('Sana`a City','Hodeidah City','Aden');
legend box off;

% Incidence

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt


load('Yemen_Gov_Incidence.mat');
ID=IData';
load('Yemen_District_Incidence.mat')
startDateofSim = datenum('10-03-2016');% Start date
XTL=datestr([startDateofSim+7.*([0:5:152])],'mm/dd/yy');

plot([1:153],ID(9,:),'k','LineWidth',2); hold on;
plot([31:153],IData(:,22),'color',hex2rgb('#CB0000'),'LineWidth',2)
plot([1:153],ID(2,:),'color',hex2rgb('#F5BE41'),'LineWidth',2)

xtickangle(45);
box off;
xlim([1 153]);
xlabel('Date','Fontsize',18);
ylabel('Suspected cholera cases','Fontsize',18);
legend('Sana`a City','Hodeidah City','Aden');
legend box off;

set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16,'XTickLabel',XTL,'XTick',([1:5:153]));

% Map
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt


S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm0_govyem_mola_20181102.shp']); % Shape file for Yemen
mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2,'Facealpha',1); hold on
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
mapshow(S(9),'FaceColor','k','Edgecolor',[0 0 0],'LineWidth',2,'Facealpha',0.35); hold on
mapshow(S(2),'FaceColor',hex2rgb('#F5BE41'),'Edgecolor',hex2rgb('#F5BE41'),'LineWidth',2,'Facealpha',0.35); hold on
SD = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen
mapshow(SD([29 31 71]),'FaceColor',hex2rgb('#CB0000'),'Edgecolor',hex2rgb('#CB0000'),'LineWidth',2,'Facealpha',0.35); hold on
axis off
mx=44.11;
MXx=44.4;
my=15.29;
Mxy=15.62;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',[0 0 0],'LineWidth',2);

mx=42.77;
MXx=43.08;
my=14.67;
Mxy=14.94;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',hex2rgb('#CB0000'),'LineWidth',2);


mx=43.34;
MXx=45.13;
my=12.57;
Mxy=12.97;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',hex2rgb('#F5BE41'),'LineWidth',2);
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm0_govyem_mola_20181102.shp']); % Shape file for Yemen

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

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2,'Facealpha',1); hold on
mapshow(SD([29 31 71]),'FaceColor',hex2rgb('#CB0000'),'Edgecolor',hex2rgb('#CB0000'),'LineWidth',2,'Facealpha',0.35); hold on

mx=42.77;
MXx=43.08;
my=14.67;
Mxy=14.94;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',hex2rgb('#CB0000'),'LineWidth',3);

axis off
xlim([42.77 43.08]);
ylim([14.67 14.94]);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2,'Facealpha',1); hold on
mapshow(SD(fS),'FaceColor','k','Edgecolor','k','LineWidth',2,'Facealpha',0.35); hold on

mx=44.11;
MXx=44.4;
my=15.29;
Mxy=15.62;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',[0 0 0],'LineWidth',3);
axis off

xlim([44.11 44.4]);
ylim([15.29 15.62]);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0808,0.18,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',2,'Facealpha',1); hold on
mapshow(SD(fA),'FaceColor',hex2rgb('#F5BE41'),'Edgecolor',hex2rgb('#F5BE41'),'LineWidth',2,'Facealpha',0.35); hold on


mx=43.34;
MXx=45.13;
my=12.57;
Mxy=12.97;
plot(linspace(mx,MXx,101),my.*ones(101,1),linspace(mx,MXx,101),Mxy.*ones(101,1),mx*ones(101,1),linspace(my,Mxy,101),MXx*ones(101,1),linspace(my,Mxy,101),'color',hex2rgb('#F5BE41'),'LineWidth',3);

xlim([43.34 45.13]);
ylim([12.57 12.97]);
axis off