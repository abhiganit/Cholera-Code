%% Governorate effectiveness of vaccinatino estimated by the model
close all;
clear;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,VT1,VT2,GNZI,GV,maxtau] = LoadYemenData; % Load the data used to construct the figure
PDS=0.85; % Specify the size of the trainign set used
atest=0; % The selction criteria

load('Fit-Vaccination-PercentData=80.mat');
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,DBE,DAE,K,n,rl,rh,mln,a,KV,dV]=RetParameterPS(par(1,:),XU,CF(:,1));

dV1=ImpactAttack(VT1(GNZI,1:length(WI(1,:))),0,dV,2,maxtau); % Two week delay until acquire immunity
dV2=ImpactAttack(VT2(GNZI,1:length(WI(1,:))),0,dV,2,maxtau);  % Two week delay until acquire immunity
EOVC=EffectOCV(dV1,KV(1),dV2,KV(2));
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

fv=find(GV(GNZI)==1);
S=S(GNZI);
startDateofSim = datenum('10-03-2016');% Start date
dW=5;
NW=length(EOVC(1,:))+maxtau;
XTL=datestr([startDateofSim+7.*[0:dW:(NW-1)]],'mm/dd/yy');


figure('units','normalized','outerposition',[0 0 1 1]);

for ii=1:sum(GV)
   
   subplot('Position',[0.045+0.5.*(rem(ii+1,2)),0.73-0.28.*(ceil(ii/2)-1),0.45,0.25]);
   bar([(1+maxtau):NW], EOVC(fv(ii),:),'k','LineStyle','none');
   box off;
   ylim([0 0.25]);
   xlim([0.5 NW+0.5]);   
   text(1, 0.24, S(fv(ii)).ADM1_EN,'Fontsize',16);
   if((ceil(ii/2)-1)==1)
    ylabel('Effectiveness of OCV campaign','Fontsize',18);
   end
   if((ceil(ii/2)-1)==2)    
        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:(NW)],'XTicklabel',XTL,'YTick',[0:0.05:0.25],'Fontsize',16,'XMinortick','on','YMinortick','on');
        xtickangle(45);
        xlabel('Week of report','Fontsize',18);
   else
        set(gca,'LineWidth',2,'tickdir','out','XTick',[1:dW:(NW)],'XTicklabel','','YTick',[0:0.05:0.25],'Fontsize',16,'XMinortick','on','YMinortick','on');
   end
end
