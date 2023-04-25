
%% Inlcudes the effwects of conflict and shellngs on the diesel prices
% Read table of past fitsclose all;
close all;
clear;
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,PopS,CI] = LoadYemenData;
CI=10000.*CI./PopS;
[~,~,CCR,~,GNZI,RC,WWRC,MI] = Contribution_to_Cholera_Incidence;
%% Plot the data
RC=RC(GNZI);

ColorM=[[152,78,163]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        [247,129,191]./255;];   % Temprature
    
WW=zeros(length(GNZI),8);
WW2=zeros(length(GNZI),7);
for mm=1:7
    temp=((CCR{mm}))./(MI);
    temp3=((CCR{mm}))./(CCR{1}+CCR{2}+CCR{3}+CCR{4}+CCR{5}+CCR{6}+CCR{7});
    for jj=1:21
        temp2=temp(jj,MI(jj,:)>0);
        WW(jj,mm) = mean(temp2);
        temp4=temp3(jj,MI(jj,:)>0);
        WW2(jj,mm) = mean(temp4);
    end
end  

mm=8;
WW(:,mm)=1-sum(WW(:,1:7),2);    
    
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

XGL={S(GNZI).ADM1_EN};
WWT=[WW(:,1:7) sum(WW(:,1:3),2)];
[WS, indexs]=sortrows(WWT,8);
WS=WS(:,1:7);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.117121848739496,0.109422492401216,0.347689075630252,0.8]);
b=barh([1:21],WS,'stacked','LineStyle','none'); 


for ii=1:7
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
               scatter(dx1(dx1<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx1(dx1<sum(WS(mm,1:7))))),5,'k','filled');
            else
                scatter(dx2(dx2<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx2(dx2<sum(WS(mm,1:7))))),5,'k','filled');
            end
        end
    end
end
xlabel('Average relative contribution per week','Fontsize',18);
h=ylabel('Governorate','Fontsize',18);
text(-0.14637462142461,23.263291139240508,'A','Fontsize',32,'Fontweight','bold');
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
xlim([10^0.5 10^4])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:21],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence per 10,000','Fontsize',18,'Color',[0.5 0.5 0.5]);

[r,p]=corr(log(CI(GNZI(indexs),end)),sum(WS(:,1:3),2),'Type','Kendall');

fprintf('Kendall Correlation between average conflict and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);


[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,4),'Type','Kendall');

fprintf('Kendall Correlation between average diesel and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);


[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,5),'Type','Kendall');

fprintf('Kendall Correlation between average wheat and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);


[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,6),'Type','Kendall');

fprintf('Kendall Correlation between average rainfall and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);

%% District
clear;

[CCR,~,~,RC,MI] = Contribution_to_Cholera_Incidence_District();

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
CI=10000.*CI./PopS;
GNZI=[4:13];
RC=RC(GNZI);

ColorM=[[152,78,163]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        [247,129,191]./255;];   % Temprature
    
WW=zeros(length(GNZI),8);
WW2=zeros(length(GNZI),7);
for mm=1:7
    temp=((CCR{mm}))./(MI);
    temp3=((CCR{mm}))./(CCR{1}+CCR{2}+CCR{3}+CCR{4}+CCR{5}+CCR{6}+CCR{7});
    for jj=1:21
        temp2=temp(jj,MI(jj,:)>0);
        WW(jj,mm) = mean(temp2);
        temp4=temp3(jj,MI(jj,:)>0);
        WW2(jj,mm) = mean(temp4);
    end
end  

mm=8;
WW(:,mm)=1-sum(WW(:,1:7),2);    
    
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
WWT=[WW(:,1:7) sum(WW(:,1:3),2)];
WWT=WWT(GNZI,:);
[WS, indexs]=sortrows(WWT,8);
WS=WS(:,1:7);


subplot('Position',[0.6,0.59,0.347689075630252,0.32]);

b=barh([1:length(GNZI)],WS,'stacked','LineStyle','none'); 
for ii=1:7
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
               scatter(dx1(dx1<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx1(dx1<sum(WS(mm,1:7))))),5,'k','filled');
            else
                scatter(dx2(dx2<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx2(dx2<sum(WS(mm,1:7))))),5,'k','filled');
            end
        end
    end
end
xlabel('Average relative contribution per week','Fontsize',18);
h=ylabel('Districts in Amanat Al Asimah','Fontsize',18);
text(-0.1632,12.592405063291142,'B','Fontsize',32,'Fontweight','bold');
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
xlim([10^0.5 10^4])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:length(GNZI)],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence per 10,000','Fontsize',18,'Color',[0.5 0.5 0.5]);


[r,p]=corr(log(CI(GNZI(indexs),end)),sum(WS(:,2:3),2),'Type','Kendall');


fprintf('(Amanat Al Asimah) Kendall Correlation between average conflict and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);



[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,4),'Type','Kendall');

fprintf('Kendall Correlation between average diesel and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);


[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,6),'Type','Kendall');

fprintf('Kendall Correlation between average rainfall and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);

%% District (Aden)
clear;



[CCR,~,~,RC,MI] = Contribution_to_Cholera_Incidence_District();

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
CI=10000.*CI./PopS;
GNZI=[14:21];
RC=RC(GNZI);

ColorM=[[152,78,163]./255; % Targeted attacks
        hex2rgb('#DE7A22'); %Conflict
        hex2rgb('#4C3F54'); %Shellings
        [153,52,4]./255; %Deisel
        hex2rgb('#FAAF08'); %Wheat
        [5,112,176]./255; %Rainfall
        [247,129,191]./255;];   % Temprature
WW=zeros(length(GNZI),8);
WW2=zeros(length(GNZI),7);
for mm=1:7
    temp=((CCR{mm}))./(MI);
    temp3=((CCR{mm}))./(CCR{1}+CCR{2}+CCR{3}+CCR{4}+CCR{5}+CCR{6}+CCR{7});
    for jj=1:21
        temp2=temp(jj,MI(jj,:)>0);
        WW(jj,mm) = mean(temp2);
        temp4=temp3(jj,MI(jj,:)>0);
        WW2(jj,mm) = mean(temp4);
    end
end  

mm=8;
WW(:,mm)=1-sum(WW(:,1:7),2);    
    
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
WWT=[WW(:,1:7) sum(WW(:,1:3),2)];
WWT=WWT(GNZI,:);
[WS, indexs]=sortrows(WWT,8);
WS=WS(:,1:7);



subplot('Position',[0.6,0.109422492401216,0.347689075630252,0.32]);
labels = {'Target attacks','Weekly conflict','Shellings/attacks','Diesel','Wheat','Rainfall','Temprature'};
b=barh([1:length(GNZI)],WS,'stacked','LineStyle','none'); 
for ii=1:7
   b(ii).FaceColor=ColorM(ii,:); 
   text(0.355,8-0.85.*(ii-1), labels{ii},'Fontsize',20,'Color',ColorM(ii,:));
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
               scatter(dx1(dx1<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx1(dx1<sum(WS(mm,1:7))))),5,'k','filled');
            else
                scatter(dx2(dx2<sum(WS(mm,1:7))),(0.15.*(ii-3)+ mm).*ones(size(dx2(dx2<sum(WS(mm,1:7))))),5,'k','filled');
            end
        end
    end
end
xlabel('Average relative contribution per week','Fontsize',18);
h=ylabel('Districts in Aden','Fontsize',18);
text(-0.1632,10.111392405063292,'C','Fontsize',32,'Fontweight','bold');
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
xlim([10^0.5 10^4])
set(ax2,'LineWidth',2,'tickdir','out','XScale','log','YTick',[1:length(GNZI)],'YTickLabel',{},'Fontsize',16);
ax2.XColor=[0.5 0.5 0.5];
xlabel('Cumulative incidence per 10,000','Fontsize',18,'Color',[0.5 0.5 0.5]);

[r,p]=corr(log(CI(GNZI(indexs),end)),sum(WS(:,2:3),2),'Type','Kendall');

fprintf('(Aden) Kendall Correlation between average conflict and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);



[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,4),'Type','Kendall');

fprintf('Kendall Correlation between average diesel and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);


[r,p]=corr(log(CI(GNZI(indexs),end)),WS(:,6),'Type','Kendall');

fprintf('Kendall Correlation between average rainfall and cumulative incidence: r= %4.3f and p=%3.2E \n',[r p]);

print(gcf,'Figure4','-depsc','-r600');
print(gcf,'Figure4','-dpng','-r600');