% Read table of past fitsclose all;
%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
NatIData=sum(IData,2);
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData'; % Transpose the data set such that the number of areas is the row
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
load('Attack_Yemen_Time_Location_Project_Forward.mat'); % Load the attacks in the area for the projection
tA=GLevelConflict(ProA,S,153);% % Send attack data (time, latitude, longatude) and shapefile of area wanted to 
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas %Load the rainfall for the projection
Rtv=Rtv(:,1:length(tA(1,:))); % truncate the railfall data to the time period spceified

% Load population density
load('PopulationDensity_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
P=log(P)';
% Load Rebel Control
load('RebelControl_Yemen.mat');
RC=RC';
% Load Health facility density
HS = shaperead([ pwd '\ShapeFile\healthsites.shp']); % Shape file for Yemen
load('PopulationSize_Yemen.mat');
H= GLevelHealthSites(HS,S);
H=10000.*H./AP';

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference

    AF=2;
    %% Forward selection
    XU=[0     1     1     1     1     1     1     0     1     0     0];
    par=[-15.9900000000000,-1.55640624891769,-0.903991753304777,-0.320435872238888,-0.930602200863455,-0.989623309352588,-2.24924929137942,-15.9900000000000,-2.27352940689430,-15.9900000000000,-15.9900000000000,0.740000000000000,0.490000000000000,0.490000000000000,0.740000000000000,0.240000000000000,0.990000000000000,0.656666666666667,0.990000000000000,-0.521594815251887,-0.0687379602147109,-8.10153681596401e-06,-1.34100092740471,0.686051158054122,-6.22022461207972,-10.5671918136867,-3.74071938385970];
    %% Complex model
    %par=[1.01805008875442,-1.58238552075403,-0.914828190362037,-0.317750487287607,-0.982399551608217,-0.996984195737628,-2.43063722171153,0.532150577454823,-2.37062441770914,0.842555841152295,-1.86758315002780,0.490000000000000,0.490000000000000,0.490000000000000,0.740000000000000,0.990000000000000,0.990000000000000,0.562916666666667,0.656666666666667,-14.6758742238635,-0.0776527356016243,-0.000198969747455813,-1.77103822419636,-13.9126220282602,-13.2014775342388,-2.16240217346506e-06,-14.0998949084456];
    %XU=ones(1,11);
    %% Adjust ascpects of functions and data for the fitting

    %% Adjust ascpects of functions and data for the fitting

    maxtau=4; % The maximum lag allowed for the model

    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
    
    %% Run the projection
    
    %% Run the logistic model with the data

    [Yt,X]= LogisticModel(beta,WI(GNZI,:),tA(GNZI,1:length(WI(1,:))),DB,DA,DBE,DAE,Ctv(GNZI,1:length(WI(1,:))),K,n,Rtv(GNZI,1:length(WI(1,:))),RIF,rl,RF,rh,tau,maxtau,CF,P(GNZI),RC(GNZI),H(GNZI));

    CInc=(beta(2).*squeeze(X(2,:,:))+beta(3).*squeeze(X(3,:,:))+beta(4).*squeeze(X(4,:,:))+beta(9).*squeeze(X(9,:,:)));
    CRain=(beta(7).*squeeze(X(7,:,:))+beta(8).*squeeze(X(8,:,:)));
    CConflict=(beta(5).*squeeze(X(5,:,:))+beta(6).*squeeze(X(6,:,:))+beta(10).*squeeze(X(10,:,:))+beta(11).*squeeze(X(11,:,:)));
    
    Mt=sum(Yt);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.0708,0.163120567375887,0.897162184873949,0.793313069908819]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

%Tot=repmat((sum(CInc,1)+sum(CRain,1)+sum(CConflict,1))./(Mt),3,1);
h=bar([1:149]+maxtau,([sum(CInc,1);sum(CConflict,1);sum(CRain,1)])','stacked','LineStyle','none');
h(1).FaceColor = 'flat';
h(1).CData = [0.4 0.4 0.4];
h(3).FaceColor = 'flat';
h(3).CData = [0 0.6 1];
h(2).FaceColor = 'flat';
h(2).CData = [0.9 0 0];
% The size to separate the weeks in the x-label

hold on;
hmf=plot([1:149]+maxtau,NatIData(1+maxtau:end),'k-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
dW=4;
NW=length(NatIData);
startDateofSim = datenum('10-03-2016');% The week of our first data point (October 3, 2016)
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:5000:55000],'Yminortick','on','Fontsize',16,'XTickLabel',XTL);

    % Sets the y-axis to not have 10^n
    ax=gca; % finds the current axis
    ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
box off;
xtickangle(45);

xlabel('Date','Fontsize',18);
yh=ylabel('Suspected cholera cases','Fontsize',18);
ylim([0 55005]);
xlim([1 153.5]);
IW=[1 21; 22 74; 75 116; 117 126];
WN=struct('N',{'First wave','Second wave','Third wave','Fourth wave'});
for ii=1:3
    plot((maxtau+mean([IW(ii,2) IW(ii+1,1)])).*ones(1001,1),linspace(0,55005,1001),'k-.','LineWidth',2);
    text((mean(IW(ii,:))),56020,WN(ii).N,'Fontsize',18);
end
text((mean(IW(4,:))),56020,WN(4).N,'Fontsize',18);

legend([h hmf],{'Incidence','Conflict','Rainfall','Data'},'Fontsize',18,'location','northwest');

legend boxoff;
text(yh.Extent(1),55528,'A','Fontsize',32,'FontWeight','bold');
print(gcf,[pwd '\Figures\Figure1A.png'],'-dpng','-r600');
%% Waves
figure('units','normalized','outerposition',[0 0 1 1]);
subplot('Position',[0.075105042016807,0.163120567375887,0.259453781512605,0.8217]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

Tot1=sum(sum(CInc(:,IW(1,1):IW(1,2)),1)+sum(CRain(:,IW(1,1):IW(1,2)),1)+sum(CConflict(:,IW(1,1):IW(1,2)),1));
Tot2=sum(sum(CInc(:,IW(2,1):IW(2,2)),1)+sum(CRain(:,IW(2,1):IW(2,2)),1)+sum(CConflict(:,IW(2,1):IW(2,2)),1));
Tot3=sum(sum(CInc(:,IW(3,1):IW(3,2)),1)+sum(CRain(:,IW(3,1):IW(3,2)),1)+sum(CConflict(:,IW(3,1):IW(3,2)),1));
Tot4=sum(sum(CInc(:,IW(4,1):IW(4,2)),1)+sum(CRain(:,IW(4,1):IW(4,2)),1)+sum(CConflict(:,IW(4,1):IW(4,2)),1));
h=barh([1:4],flip([sum(sum(CConflict(:,IW(1,1):IW(1,2)),1))./Tot1 sum(sum(CRain(:,IW(1,1):IW(1,2)),1))./Tot1;sum(sum(CConflict(:,IW(2,1):IW(2,2)),1))./Tot2 sum(sum(CRain(:,IW(2,1):IW(2,2)),1))./Tot2;sum(sum(CConflict(:,IW(3,1):IW(3,2)),1))./Tot3 sum(sum(CRain(:,IW(3,1):IW(3,2)),1))./Tot3;sum(sum(CConflict(:,IW(4,1):IW(4,2)),1))./Tot4 sum(sum(CRain(:,IW(4,1):IW(4,2)),1))./Tot4]),'stacked','LineStyle','none');

h(2).FaceColor = 'flat';
h(2).CData = [0 0.6 1];
h(1).FaceColor = 'flat';
h(1).CData = [0.9 0 0];
legend({'Conflict','Rainfall'},'Fontsize',18);
legend boxoff;
N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});

% changing the aspects of the axis for the the current figure 
set(gca,'Tickdir','out','LineWidth',2,'YTick',[1:4],'XTick',[0:0.1:1],'Xminortick','on','Fontsize',16,'YTickLabel',{'Fourth wave','Third wave','Second wave','First wave'});

box off;
xtickangle(45);

xlabel({'Estimated contribtuion to suscpeted','cholera cases'},'Fontsize',18);
xlim([0 1]);
ylim([0.5 4.5]);
text(-0.290688259109312,4.443896424167694,'B','Fontsize',32,'FontWeight','bold');

%% Gov level
%figure('units','normalized','outerposition',[0 0 1 1]);
 % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt

% Tot=repmat(sum(CInc,2)+sum(CRain,2)+sum(CConflict,2),1,2);
% h=bar([1:20],([sum(CConflict,2) sum(CRain,2) ]./Tot),'stacked','LineStyle','none');
% 
% h(2).FaceColor = 'flat';
% h(2).CData = [0 0.6 1];
% h(1).FaceColor = 'flat';
% h(1).CData = [0.9 0 0];
% legend({'Conflict','Rainfall'},'Fontsize',18);
% legend boxoff;
% N=struct('G',{'Abyan';'Aden';'Al Bayda';'Al Dhale''e';'Al Hudaydah';'Al Jawf';'Al Maharah';'Al Mahwit';'Amanat Al Asimah';'Amran';'Dhamar';'Hadramaut';'Hajjah';'Ibb';'Lahj';'Marib';'Raymah';'Sa''ada';'Sana''a';'Shabwah';'Socotra';'Taizz'});
% 
% % changing the aspects of the axis for the the current figure 
% set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',[0:0.2:1],'Yminortick','on','Fontsize',16,'XTickLabel',{N(GNZI).G});
% 
% box off;
% xtickangle(45);
% 
% xlabel('Governorate','Fontsize',18);
% ylabel({'Estimated contribtuion to suscpeted','cholera cases'},'Fontsize',18);
% ylim([0 1.02]);
% xlim([0.5 20.5])

%% Govnorate level various waves and the impact of rainfall
% 
Tot1=(sum(CInc(:,IW(1,1):IW(1,2)),2)+sum(CRain(:,IW(1,1):IW(1,2)),2)+sum(CConflict(:,IW(1,1):IW(1,2)),2));
Tot2=(sum(CInc(:,IW(2,1):IW(2,2)),2)+sum(CRain(:,IW(2,1):IW(2,2)),2)+sum(CConflict(:,IW(2,1):IW(2,2)),2));
Tot3=(sum(CInc(:,IW(3,1):IW(3,2)),2)+sum(CRain(:,IW(3,1):IW(3,2)),2)+sum(CConflict(:,IW(3,1):IW(3,2)),2));
Tot4=(sum(CInc(:,IW(4,1):IW(4,2)),2)+sum(CRain(:,IW(4,1):IW(4,2)),2)+sum(CConflict(:,IW(4,1):IW(4,2)),2));
% 
TT=[Tot1 Tot2 Tot3 Tot4];
R=[sum(CRain(:,IW(1,1):IW(1,2)),2)./Tot1 sum(CRain(:,IW(2,1):IW(2,2)),2)./Tot2 sum(CRain(:,IW(3,1):IW(3,2)),2)./Tot3 sum(CRain(:,IW(4,1):IW(4,2)),2)./Tot4];
CM=[sum(CConflict(:,IW(1,1):IW(1,2)),2)./Tot1 sum(CConflict(:,IW(2,1):IW(2,2)),2)./Tot2 sum(CConflict(:,IW(3,1):IW(3,2)),2)./Tot3 sum(CConflict(:,IW(4,1):IW(4,2)),2)./Tot4;];
for jj=1:4
    subplot('Position',[0.38,0.163120567375887+0.21*(jj-1),0.57,0.189]);
    yyaxis left
    h=bar([1:20],[CM(:,(5-jj)) R(:,(5-jj))],'stacked','LineStyle','none');
    h(2).FaceColor = 'flat';
    h(2).CData = [0 0.6 1];
    h(1).FaceColor = 'flat';
    h(1).CData = [0.9 0 0];
    box off;
    xlim([0.5 20.5]);
    ylim([0 1]);
    if(jj==1)        
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',[0:0.2:1],'Yminortick','on','Fontsize',16,'XTickLabel',{N(GNZI).G},'ycolor','k');
        xlabel('Govnorate','Fontsize',18)   
        xtickangle(45);
    else
        if(jj==2)            
            yh=ylabel({'Estimated contribtuion to suscpeted cholera cases'},'Fontsize',18);
            yh.Position=[-0.416870970915301,1.002674273628602,-1];
        end
        if(jj==4)
            text(-1,0.9591,'C','Fontsize',32','FontWeight','bold');
        end
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',[0:0.2:1],'Yminortick','on','Fontsize',16,'XTickLabel',{''},'ycolor','k');
    end
    yyaxis right
    semilogy([1:20],TT(:,5-jj),'-o','color',[0 0 0],'LineWidth',2,'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
    if(jj==1)        
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',10.^[0:6],'Fontsize',16,'XTickLabel',{N(GNZI).G},'ycolor',[0 0 0]);
        xlabel('Govnorate','Fontsize',18)   
        xtickangle(45);
    else        
        if(jj==2)  
            yh=ylabel({'Estimated suscpeted cholera cases'},'Fontsize',18);
            yh.Rotation=270;
            yh.Position=[22.13655849975374,3383877.445489365,-1];
        end
        set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:20],'YTick',10.^[0:6],'Fontsize',16,'XTickLabel',{''},'ycolor',[0 0 0]);
    end
    box off;
    
    xlim([0.5 20.5]);
    ylim([1,10^6]);
end

print(gcf,[pwd '\Figures\Figure1B.png'],'-dpng','-r600');


% Sm=cell(length(S),1); % cell array for the names of areas in shape file
% for ii=1:length(Sm)
%     Sm{ii}=S(ii).ADM1_EN; % record the names of the shape file
% end
% 
% PLab=struct('N',{'i)','ii)','iii)','iv)'});
% for ww=1:4
%     f1=figure('units','normalized','outerposition',[0 0 1 1]);
% 
%     for ii=1:length(GNZI)    
%         ff=find(contains(Sm,N(GNZI(ii)).G)); % see if the name is in the shape file based on the ordering OG
%         mapshow(S(ff),'FaceColor',[0 0 0.6],'Edgecolor',[0 0 0],'LineWidth',1,'Facealpha',R(ii,ww));hold on % shape the area in with color C(ii,:)
%     end
% 
%     mapshow(S(12),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
% 
%     mapshow(S(21),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
%     xl=xlim; % saves the xlim to xl
%     yl=ylim;% saves the ylim to yl
% 
%     % Produces a cross hatch for polygon 12
%     CrossHatchMap(min(S(12).X),max(S(12).X),min(S(12).Y),max(S(12).Y),0.01,[0 0 0],S(12))
%     % Produces a cross hatch for polygon 21
%     CrossHatchMap(min(S(21).X),max(S(21).X),min(S(21).Y),max(S(21).Y),0.025,[0 0 0],S(21))
%     xlim(xl); % set xlimits to xl
%     ylim(yl); % set ylimits to yl
%     set(gca,'visible','off') % remove the axis from the map 
%     if(ww==1)
%         text(min(xl),max(yl),'B','Fontsize',32,'Fontweight','bold');
%     end
%     text(41.8,max(yl),PLab(ww).N,'Fontsize',28);    
%     print(gcf,[pwd '\Figures\Figure1B_' num2str(ww) '.png'],'-dpng','-r600');
% end
% 
% 
% %% Govnorate level various waves and the impact of conflict
% 
% 
% 
% for ww=1:4
%     f1=figure('units','normalized','outerposition',[0 0 1 1]);
% 
%     for ii=1:length(GNZI)    
%         ff=find(contains(Sm,N(GNZI(ii)).G)); % see if the name is in the shape file based on the ordering OG
%         mapshow(S(ff),'FaceColor',[0.6 0 0],'Edgecolor',[0 0 0],'LineWidth',1,'Facealpha',CM(ii,ww));hold on % shape the area in with color C(ii,:)
%     end
% 
%     mapshow(S(12),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
% 
%     mapshow(S(21),'FaceColor',[1 1 1],'Edgecolor',[0 0 0],'LineWidth',1)
%     xl=xlim; % saves the xlim to xl
%     yl=ylim;% saves the ylim to yl
% 
%     % Produces a cross hatch for polygon 12
%     CrossHatchMap(min(S(12).X),max(S(12).X),min(S(12).Y),max(S(12).Y),0.01,[0 0 0],S(12))
%     % Produces a cross hatch for polygon 21
%     CrossHatchMap(min(S(21).X),max(S(21).X),min(S(21).Y),max(S(21).Y),0.025,[0 0 0],S(21))
%     xlim(xl); % set xlimits to xl
%     ylim(yl); % set ylimits to yl
%     set(gca,'visible','off') % remove the axis from the map 
%     if(ww==1)
%         text(min(xl),max(yl),'C','Fontsize',32,'Fontweight','bold');
%     end
%     text(41.8,max(yl),PLab(ww).N,'Fontsize',28);    
%     print(gcf,[pwd '\Figures\Figure1C_' num2str(ww) '.png'],'-dpng','-r600');
% end
