%% Compare District level incidence of the various models
%close all;
clear;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,~,GNZI,maxtau] = LoadYemenDistrictData;
PDS=0.8;
atest=0;
%% Forward selection

load(['ForwardSelectionNoRain-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
XU=XUv(end,:);
par=parv(end,:);

% Evaluate the number of paramters that are being used in the estimation 
[k,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,CF,RIF,RF]=RetParameterPS(par,XU);
%% Run the projection

%% Run the logistic model with the data

[Yt,~]= LogisticModel(beta,WI,tA,DB,DA,DBE,DAE,Ctv,K,n,Rtv,RIF,rl,RF,rh,tau,maxtau,CF,P,RC,H,WPIN,Mt,[]);
load('PopulationSize_DistrictYemen.mat');
NW2016=ceil((datenum('12-31-2016')-datenum('10-03-2016'))./7); % Number of weeks to rpelicate populatino density for 2016
NW2019=153-52-52-NW2016; % Numebr of weeks to reproduce population density for 2019
PopS=[ repmat(AP(:,1),1,NW2016) repmat(AP(:,2),1,52)  repmat(AP(:,3),1,52)  repmat(AP(:,4),1,NW2019)]; 
PopS=PopS(:,31:end); % Incidence data for the districts starts at may 1 2017
Yt=(Yt./(10000)).*PopS(:,maxtau+1:end);

load('Yemen_District_Incidence.mat'); % Incidence data
WI=IData'; % Transpose the data set such that the number of areas is the row


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

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

SName=cell(length(SD)+3,1);

for ii=1:length(SD)
   SName{ii}=SD(ii).ADM2_EN; 
end

SName{length(SD)+1}=string({'Hodeidah City'});
SName{length(SD)+2}=S(9).ADM1_EN; 
SName{length(SD)+3}=S(2).ADM1_EN; 
for ii=1:ceil(length(SName)/4)    
    figure('units','normalized','outerposition',[0 0 1 1]); % opens a new figure
    for jj=1:4
        if((jj+4*(ii-1))<=length(WI(:,1)))
            subplot('Position',[0.0708+rem(jj+1,2).*0.97/2,(0.9/2+0.163120567375887)-floor(jj/3).*0.9/2,0.897162184873949/2*0.95,0.793313069908819/2*0.82]); % Creates a sub-panel to plot the figure in the x position is 0.0708 the y position is 0.163120567375887, 0.897162184873949 is the width and 0.793313069908819 is the heigt
            D=WI((jj+4*(ii-1)),:); % Pull weekly incidence data for governonate

            MI=floor(log10(max(Yt(jj+4*(ii-1),:)))); % The max magnitude we want for the figure
            FF=ceil(max(Yt(jj+4*(ii-1),:))./10^MI);
            NW=length(D); % Number of weeks
            scatter([1:NW],D,40,'k','filled'); hold on; % Plot data
            plot([(1+(length(D)-length(Yt(jj+4*(ii-1),:)))):NW],Yt(jj+4*(ii-1),:),'b','LineWidth',2);  % Plot model predictions  

            box off; % removes the outside box on the figure
            xlim([1 NW]); % sets the x-limits of our x -axis
            ylim([0 FF*10^MI]); %sets the y-limts of the y-axis
            % The size to separate the weeks in the x-label
            dW=7;
            % Set the X-tick labels to be Dates rather than numbers
            startDateofSim = datenum('05-01-2017');% The week of our first data point (October 3, 2016)
            XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');
            % changing the aspects of the axis for the the current figure 
            if(jj>2)
                set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:10^MI:FF*10^MI],'Fontsize',16,'XTickLabel',XTL);
            % Xlable of the figure
            xlabel('Week of report','Fontsize',24);
            else
                set(gca,'Tickdir','out','LineWidth',2,'XTick',[1:dW:NW],'YTick',[0:10^MI:FF*10^MI],'Fontsize',16,'XTickLabel',{});
            end
            % Rotates the xticklabel 
            xtickangle(45) 
            % Sets the y-axis to not have 10^n
            ax=gca; % finds the current axis
            ax.YAxis.Exponent = 0; % Sets the y-axis to not have 10^n
            title(SName(jj+4*(ii-1)),'Fontsize',28);
            if(rem(jj,2)==1)
            % ylable of the figure
            ylabel('Number of reported cases','Fontsize',24);
            end
        end
    end
end