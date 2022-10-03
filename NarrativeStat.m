close all;
clear;
clc;
[WID,~,tAd,RtvD,~,PD,RCD,HD,WPINmD,FPINmD,DieseltD,WheattD,V1D,V2D,GNZID,GVD,~,~,CID] = LoadYemenDistrictData; % Load the data used to construct the figure
[WI,~,tA,Rtv,~,P,RC,H,WPINm,FPINm,Dieselt,Wheatt,V1,V2,GNZI,GV,~,~,CI] = LoadYemenData; % Load the data used to construct the figure
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

% Record the names for the IDP calculation
Sm=cell(length(S),1); % allocate space
for ii=1:length(Sm)
    Sm{ii}=S(ii).ADM1_EN; % record name
end

load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
CtG=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
MtG=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize



load('Conflict_Yemen_Time_Location_Prior_to_Cholera_Epidemic.mat')% Load the conflict in the area for the projection

BeforeCholeraStart=abs(min(ProCBefore(:,1)))+1;
ProCBefore(:,1)=ProCBefore(:,1)+abs(min(ProCBefore(:,1)))+1;


CtGBefore=GLevelConflict(ProCBefore,S,BeforeCholeraStart);
load('Yemen_Air_Shelling_Before_Cholera.mat');

BeforeCholeraStart=abs(min(YAStBefore(:,1)))+1;
YAStBefore(:,1)=YAStBefore(:,1)+abs(min(YAStBefore(:,1)))+1;
MtGBefore=GLevelConflict(YAStBefore,S,BeforeCholeraStart); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize


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

load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
CtD=GLevelConflict(ProC,SD,153);

temp=GLevelConflict(ProC,S([9 2]),153);
CtD=[CtD; sum(CtD(1:3,:),1);temp];

load('Yemen_Air_Shelling.mat');
MDt=GLevelConflict(YASt,SD,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
temp=GLevelConflict(YASt,S([9 2]),153);
MDt=[MDt; sum(MDt(1:3,:),1);temp];


load('Conflict_Yemen_Time_Location_Prior_to_Cholera_Epidemic.mat')% Load the conflict in the area for the projection

BeforeCholeraStart=abs(min(ProCBefore(:,1)))+1;
ProCBefore(:,1)=ProCBefore(:,1)+abs(min(ProCBefore(:,1)))+1;


CtDGBefore=GLevelConflict(ProCBefore,SD,BeforeCholeraStart);

temp=GLevelConflict(ProCBefore,S([9 2]),BeforeCholeraStart);

CtDGBefore=[CtDGBefore; sum(CtDGBefore(1:3,:),1);temp];

load('Yemen_Air_Shelling_Before_Cholera.mat');

BeforeCholeraStart=abs(min(YAStBefore(:,1)))+1;
YAStBefore(:,1)=YAStBefore(:,1)+abs(min(YAStBefore(:,1)))+1;
MDtGBefore=GLevelConflict(YAStBefore,SD,BeforeCholeraStart); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

temp=GLevelConflict(YAStBefore,S([9 2]),BeforeCholeraStart);
MDtGBefore=[MDtGBefore; sum(MDtGBefore(1:3,:),1);temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Sana'a City
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

fprintf('================================================================================ \n')
fprintf([Sm{9} '\n']);
fprintf('================================================================================ \n \n')

fprintf(['Shelling and air attack in ' Sm{9} ' from Janurart 1, 2015 to Oct 2 2016 = %d \n'],sum(MtGBefore(9,:)));

WeekEnd=floor((datenum('April 9, 2017')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('October 3, 2016')-datenum('October 3, 2016'))./7)+1;

fprintf(['Shelling and air attack in ' Sm{9} ' between April 9, 2017 and October 3, 2016 = %d \n'],sum(MtG(9,WeekStart:WeekEnd)));


WeekEnd2=floor((datenum('November 20, 2016')-datenum('October 3, 2016'))./7)+1;

fprintf(['Percentage of shelling and air attack in ' Sm{9} ' between November 20, 2016 and October 3, 2016 = %3.1f %% \n'],100.*sum(MtG(9,WeekStart:WeekEnd2))/sum(MtG(9,WeekStart:WeekEnd)));



WeekEnd=floor((datenum('April 23, 2017')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('October 3, 2016')-datenum('October 3, 2016'))./7)+1;

MD=sum(MDt(4:(4+length(fS)-1),WeekStart:WeekEnd),2);
ND=cell(length(fS),1);
for ii=1:length(ND)
    ND{ii}=SD(4+(ii-1)).ADM2_EN;
end

fprintf(['Maximum shelling and air attack in ' ND{MD==max(MD)} ' between April 23, 2017 and October 3, 2016 = %d \n'],max(MD));
fprintf(['Minimum shelling and air attack in ' ND{MD==min(MD)} ' between April 23, 2017 and October 3, 2016 = %d \n'],min(MD));


%Second wave: March 27, 2017 to April 1, 2018-


WeekEnd=floor((datenum('April 1, 2018')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('March 27, 2017')-datenum('October 3, 2016'))./7)+1;

FM=MtG(9,WeekStart:WeekEnd);
fnz=find(FM>0);
fnz= fnz(2:end)- fnz(1:end-1);
fprintf(['Frequency of at least one shelling per week in ' Sm{9} ' between March 27, 2017 and April 1, 2018 = %3.2f \n'],mean(fnz));

WeekEnd=floor((datenum('April 29, 2018')-datenum('May 1, 2017'))./7)+1;
WeekStart=floor((datenum('May 1, 2017')-datenum('May 1, 2017'))./7)+1;

MA=mean(WID(4:(4+length(fS)-1),WeekStart:WeekEnd),2);

fprintf(['Average attack rate in ' ND{MA==max(MA)} ' between May 1, 2017 and April 29, 2018 = %3.2f per 10,000 \n'],max(MA));

WeekEnd=floor((datenum('April 1,2018')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('March 27, 2017')-datenum('October 3, 2016'))./7)+1;
fprintf(['Average attacks per week in ' Sm{9} ' in early stage vs after the second wave: %3.2f to %3.2f \n'],[mean(MtG(9,WeekStart:WeekEnd)) mean(MtG(9,(WeekEnd+1):end))])

fprintf(['Decrease attacks per week in ' SD(7).ADM2_EN ' in early stage vs after the second wave: %3.2f%% \n'],100.*[1-mean(MDt(7,(WeekEnd+1):end))./mean(MDt(7,WeekStart:WeekEnd))])

fprintf(['Increase attacks per week in ' SD(11).ADM2_EN ' in early stage vs after the second wave: %3.2f%% \n'],100.*[mean(MDt(11,(WeekEnd+1):end))./mean(MDt(11,WeekStart:WeekEnd))-1])


WeekEnd=floor((datenum('April 1, 2018')-datenum('May 1, 2017'))./7)+1;
WeekStart=floor((datenum('May 1, 2017')-datenum('May 1, 2017'))./7)+1;

fprintf(['Attack rate in ' SD(7).ADM2_EN ' second wave: %3.1f per 10000 \n'],sum(WID(7,(WeekStart):WeekEnd)));
fprintf(['Attack rate in ' SD(11).ADM2_EN ' second wave: %3.1f per 10000 \n'],sum(WID(11,(WeekStart):WeekEnd)));


fprintf(['Attack rate in ' SD(7).ADM2_EN ' after the second wave: %3.1f per 10000 \n'],sum(WID(7,(WeekEnd+1):end)));
fprintf(['Attack rate in ' SD(11).ADM2_EN ' after second wave: %3.1f per 10000 \n'],sum(WID(11,(WeekEnd+1):end)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Hodeidah City
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fprintf('\n================================================================================ \n')
fprintf(['Hodeidah City \n']);
fprintf('================================================================================ \n \n')
%Second wave: March 27, 2017 to April 1, 2018

fprintf(['Shelling and air attack in Hodeidah City vs ' Sm{9} ' from Janurart 1, 2015 to Oct 2 2016 = %d vs %d \n'],[sum(MDtGBefore(end-2,:)) sum(MtGBefore(9,:)) ]);

WeekEnd=floor((datenum('March 26, 2017')-datenum('October 3, 2016'))./7)+1;

fprintf('Attacks in Hodeidah City during the primary wave = %d \n',sum(MDt(end-2,1:WeekEnd)));
fprintf('Percentage of attacks in Hodeidah City during the primary wave = %3.2f \n',100.*sum(MDt(end-2,1:WeekEnd))./sum(MtG(5,1:WeekEnd)));

WeekEnd2=floor((datenum('March 26, 2017')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('January 16, 2017')-datenum('October 3, 2016'))./7)+1;

fprintf('Percentage of attacks in Hodeidah City during January 16, 2017 and  Marhc 26, 2017= %3.2f%% \n',100.*sum(MDt(end-2,WeekStart:WeekEnd2))./sum(MDt(end-2,1:WeekEnd)));


WeekEnd=floor((datenum('April 29, 2018')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('May 1, 2017')-datenum('October 3, 2016'))./7)+1;

fprintf('Attacks in Hodeidah City during May 1, 2017 and April 29, 2018 = %d \n',sum(MDt(end-2,WeekStart:WeekEnd)));

FM=MDt(end-2,WeekStart:WeekEnd);
fnz=find(FM>0);
fnz= fnz(2:end)- fnz(1:end-1);
fprintf(['Frequency of at least one shelling per week in Hodeidah City during May 1, 2017 and April 29, 2018 = %3.1f \n'],mean(fnz));

ND=cell(3,1);
for ii=1:length(ND)
    ND{ii}=SD(ii).ADM2_EN;
    fprintf(['Attacks in ' ND{ii} ' during May 1, 2017 and April 29, 2018 = %d \n'],sum(MDt(ii,WeekStart:WeekEnd)));
end

for ii=1:length(ND)
    fprintf(['Attack rate in ' ND{ii} ' during May 1, 2017 and April 29, 2018 = %3.2f per 10000 \n'],mean(WID(ii,WeekStart:WeekEnd)));
end

FM=MDt(end-2,(WeekEnd+1):end);
fnz=find(FM>0);
fnz= fnz(2:end)- fnz(1:end-1);
fprintf(['Frequency of at least one shelling per week in Hodeidah City after the second wave = %3.1f \n'],mean(fnz));


fprintf('mean number Attacks in Hodeidah City during May 1, 2017 and April 29, 2018 = %3.1f \n',mean(MDt(end-2,WeekStart:WeekEnd)));
fprintf('mean number Attacks in Hodeidah City after second wave = %3.1f \n',mean(MDt(end-2,(WeekEnd+1):end)));
figure('units','normalized','outerposition',[0.0 0.025 1 1]);
plot([1:length(MDt(end-2,:))],MDt(end-2,:),'k','LineWidth',2);
XTL=datestr(datenum('October 3,2016') +7.*[0:1:(length(MDt(end-2,:))-1)]);
set(gca,'tickdir','out','XTick',[1:1:length(MDt(end-2,:))],'XTickLabel',{XTL});
box off;
ylabel('Number of shelling attacks');
title('Hodeidah City');
xtickangle(45);
fprintf(['Peak shelling in in Hodeidah City ' XTL(MDt(end-2,:)==max(MDt(end-2,:)),:) '\n']);

figure('units','normalized','outerposition',[0.0 0.025 1 1]);
plot([1:length(WID(end-2,:))],WID(end-2,:),'r','LineWidth',2);
XTL2=datestr(datenum('May 1,2017') +7.*[0:1:(length(MDt(end-2,:))-1)]);
set(gca,'tickdir','out','XTick',[1:1:length(MDt(end-2,:))],'XTickLabel',{XTL});
box off;
ylabel('Attack rate');
title('Hodeidah City');
xtickangle(45);


for ii=1:length(ND)
    fprintf(['Attacks in ' ND{ii} ' after April 29, 2018 = %d \n'],sum(MDt(ii,(WeekEnd+1):end)));
end

for ii=1:length(ND)
    fprintf(['Average attack rate in ' ND{ii} ' after April 29, 2018 = %3.2f per 10000 \n'],mean(WID(ii,(WeekEnd+1):end)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Aden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fprintf('\n================================================================================ \n')
fprintf([Sm{2} '\n']);
fprintf('================================================================================ \n \n')

fprintf(['Shelling and air attack in ' Sm{2} ' from Janurart 1, 2015 to Oct 2 2016 = %d \n'],sum(MtGBefore(2,:)));

fprintf(['Perecentage of Shelling and air attack in ' Sm{2} ' during March 23 and August 2 for period Janurart 1, 2015 to Oct 2 2016 = %4.3f \n'],sum(MtGBefore(2,13:31))./sum(MtGBefore(2,:)));

WeekEnd=floor((datenum('March 26, 2017')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('October 3, 2016')-datenum('October 3, 2016'))./7)+1;
MD=sum(MDt((4+length(fS)):(4+length(fS)+length(fA)-1),WeekStart:WeekEnd),2);
ND=cell(length(fA),1);
for ii=1:length(ND)
    ND{ii}=SD(4+length(fS)+(ii-1)).ADM2_EN;
end

for ii=1:length(ND)
    fprintf(['Attacks in ' ND{ii} ' between Oct 3 2016 and March 26 2017 = %d \n'],MD(ii));
end


%Second wave: March 27, 2017 to April 1, 2018-


WeekEnd=floor((datenum('April 1, 2018')-datenum('October 3, 2016'))./7)+1;
WeekStart=floor((datenum('March 27, 2017')-datenum('October 3, 2016'))./7)+1;
MD=sum(MDt((4+length(fS)):(4+length(fS)+length(fA)-1),WeekStart:WeekEnd),2);

for ii=1:length(ND)
    fprintf(['Attacks in ' ND{ii} ' between March 27 2017 and April 1, 2018 = %d \n'],MD(ii));
end

dWID=mean(WID((4+length(fS)):(4+length(fS)+length(fA)-1),WeekStart:WeekEnd),2);
for ii=1:length(ND)
    if((dWID(ii)==max(dWID))||(dWID(ii)==min(dWID)))
        fprintf(['Average attack rate in ' ND{ii} ' during May 1, 2017 and April 1, 2018 = %3.2f per 10000 \n'],dWID(ii));
    end
end

fprintf(['Attacks in Aden after the second wave = %d \n'],sum(MtG(2,(WeekEnd+1):end)));

FM=MtG(2,(WeekEnd+1):end);
fnz=find(FM>0);
fnz= fnz(2:end)- fnz(1:end-1);
fprintf(['Frequency of at least one shelling per week in ' Sm{2} ' after the second wave = %3.2f \n'],mean(fnz));

MD=sum(MDt((4+length(fS)):(4+length(fS)+length(fA)-1),(WeekEnd+1):end),2);
for ii=1:length(ND)
    fprintf(['Attacks in ' ND{ii} ' after the second wave = %d \n'],MD(ii));
end

dWID=mean(WID((4+length(fS)):(4+length(fS)+length(fA)-1),(WeekEnd+1):end),2);
for ii=1:length(ND)
    if((dWID(ii)==max(dWID))||(dWID(ii)==min(dWID)))
        fprintf(['Average attack rate in ' ND{ii} ' after second wave = %3.2f per 10000 \n'],dWID(ii));
    end
end