%% Plots the location of the attacks on water related systme before and after cholera cases were reported
clc;
clear;
close all;
load('Yemen_Air_Shelling.mat');
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


%% After the outbreak reprted cases
f1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.143120567375887,0.897162184873949,0.793313069908819]);

mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on
f=find(YASt(:,1)>=1); % find the attacks that have happened before 
LLB=YASt(f,2:3);
B=unique(LLB,'rows');
B=[B zeros(length(B(:,1)),1)];

for ii=1:length(B(:,1))
    f=find(B(ii,1)==LLB(:,1));
    g=find(B(ii,2)==LLB(f,2));
    B(ii,3)=length(g);
end

% scatter(B(:,2),B(:,1),(B(:,3)),'k','filled')

for ii=1:length(B(:,1))

    t = linspace(0, 2*pi);
    r = 0.05;
    x = r*cos(t)+B(ii,1);
    y = r*sin(t)+B(ii,2);
   h=patch(y,x, 'k','Facealpha',0.35,'EdgeAlpha',0) ;
end
scatter(53.08.*ones(5,1),[16.5 16.2 15.9 15.6 15.3]+2.5,[500 100 50 25 5],'k','filled')
axis off;
text(53.25.*ones(5,1),[16.5 16.2 15.9 15.6 15.3]+2.5,{'500 Attacks','100 Attacks','50 Attacks','25 Attacks','5 Attacks'},'Fontsize',16);
load('Attack_Yemen_Time_Location_Project_Forward.mat');
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

%% After the outbreak reprted cases
f1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.143120567375887,0.897162184873949,0.793313069908819]);

mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on
f=find(ProA(:,1)>=1); % find the attacks that have happened before 
LLB=ProA(f,2:3);
B=unique(LLB,'rows');
B=[B zeros(length(B(:,1)),1)];

for ii=1:length(B(:,1))
    f=find(B(ii,1)==LLB(:,1));
    g=find(B(ii,2)==LLB(f,2));
    B(ii,3)=length(g);
end

scatter(B(:,2),B(:,1),62.5.*(B(:,3)),'r','filled')

scatter(53.08.*ones(5,1),[16.5 16.2 15.9 15.6 15.3]+2.5,62.5.*[8 6 4 2 1],'r','filled')
axis off;
text(53.25.*ones(5,1),[16.5 16.2 15.9 15.6 15.3]+2.5,{'8 Attacks','6 Attacks','4 Attacks','2 Attacks','1 Attack'},'Fontsize',16,'color','r');
