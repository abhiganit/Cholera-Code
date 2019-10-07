%% Plots the location of the attacks on water related systme before and after cholera cases were reported
clc;
clear;
close all;
load('AllAttack_Yemen_Time_Location_Project_Forward.mat');
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

%% Before the outbreak
f1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.143120567375887,0.897162184873949,0.793313069908819]);

mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on
f=find(ProCA(:,1)<1); % find the attacks that have happened before 
LLB=ProCA(f,2:3);
B=unique(LLB,'rows');
B=[B zeros(length(B(:,1)),1)];

for ii=1:length(B(:,1))
    f=find(B(ii,1)==LLB(:,1));
    g=find(B(ii,2)==LLB(f,2));
    B(ii,3)=length(g);
end

scatter(B(:,2),B(:,1),30.*B(:,3),'k','filled')
axis off;
title('Attacks on water systems before cases reported','Fontsize',32);

%% After the outbreak reprted cases
f1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.143120567375887,0.897162184873949,0.793313069908819]);

mapshow(S,'FaceColor','none','Edgecolor',[0 0 0],'LineWidth',1); hold on
f=find(ProCA(:,1)>=1); % find the attacks that have happened before 
LLB=ProCA(f,2:3);
B=unique(LLB,'rows');
B=[B zeros(length(B(:,1)),1)];

for ii=1:length(B(:,1))
    f=find(B(ii,1)==LLB(:,1));
    g=find(B(ii,2)==LLB(f,2));
    B(ii,3)=length(g);
end

scatter(B(:,2),B(:,1),30.*B(:,3),'r','filled')
axis off;
title('Attacks on water systems after cases reported','Fontsize',32);
