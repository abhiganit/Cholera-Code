close all;
clear;
FC=hex2rgb('#a6bddb');

[Y_Avg,NW,NW1,IData_VAL,maxtau,GNZI]  = Weighted_Model_Incidence('District');

%% District

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

for ii=1:(length(SD)) % Do not want to examine the last three as these are the larger areas and not the districts
    startDateofSim = datenum('5-01-2016');% Start date
    figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

    subplot('Position',[0.07,0.225,0.92,0.755]);
    dW=10;
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    
    bar([(1+maxtau):NW],Y_Avg(ii,:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
    scatter([1:NW],IData_VAL(ii,:),20,'k','filled'); 
    box off;
    xlim([0.5 NW+0.5]);
    set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTick',[0:1000:8000],'YMinortick','on');
    ylim([0 max(ylim).*1.1])
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    xtickangle(90);
    hy=ylabel('Suspected cholera cases','Fontsize',18);
%     ylim([0 100]);

text(1.5,max(ylim),['\bf{\it{' SD(ii).ADM2_EN '}}'],'Fontsize',18);
    text(1.5,max(ylim)*(1-0.06),'Model estimate','Fontsize',18,'Color',FC(1,:));
%     text(1.5,max(ylim)*(1-0.12),'Model validation','Fontsize',18,'Color',FC(2,:));
    text(1.5,max(ylim)*(1-0.12),'Reported cases','Fontsize',18);
    

    xlabel('Week reported','Fontsize',18);
    
    print(gcf,['District_Validation_' SD(ii).ADM2_EN '.png'],'-dpng','-r600');
end

ii=22;
% Al Hudaydah City
figure('units','normalized','outerposition',[0.05 0.05 0.8 0.65]);

    subplot('Position',[0.07,0.225,0.92,0.755]);
    dW=10;
    XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
    
    bar([(1+maxtau):NW],Y_Avg(ii,:),'Facecolor',FC,'LineStyle','none','Facealpha',0.6); hold on
    scatter([1:NW],IData_VAL(ii,:),20,'k','filled'); 
    box off;
    xlim([0.5 NW+0.5]);
    set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:NW],'XTickLabel',XTL,'Fontsize',16,'Xminortick','on');%,'YTick',[0:1000:8000],'YMinortick','on');
    ylim([0 max(ylim).*1.1])
    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
    xtickangle(90);
    hy=ylabel('Suspected cholera cases','Fontsize',18);
%     ylim([0 100]);

text(1.5,max(ylim),['\bf{\it{Al Hudaydah City}}'],'Fontsize',18);
    text(1.5,max(ylim)*(1-0.06),'Model estimate','Fontsize',18,'Color',FC(1,:));
%     text(1.5,max(ylim)*(1-0.12),'Model validation','Fontsize',18,'Color',FC(2,:));
    text(1.5,max(ylim)*(1-0.12),'Reported cases','Fontsize',18);
    

    xlabel('Week reported','Fontsize',18);
    
    print(gcf,['District_Validation_Hodeidah_City.png'],'-dpng','-r600');