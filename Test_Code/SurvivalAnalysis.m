%% Prints a summary of the information for the different govnorates
clear;
clc;
close all;
% 
%% Cumulative Cases

IW=[1 21; 22 74; 75 116; 117 153; 1 153];
PNS=10000;
TES=PNS.*[0.0005 0.01 0.001 0.005 0.025];
PT=zeros(5,6);
for Wv=1:5
    load('Yemen_Gov_Incidence.mat'); % Incidence data
    load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
    S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
    load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
    WI=IData; % Transpose the data set such that the number of areas is the row
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

    load('PopulationSize_Yemen.mat')
    load('Yemen_Gov_Incidence.mat')
    wc=cumsum(IData(IW(Wv,1):IW(Wv,2),:));
    PP=repmat(AP,IW(Wv,2)-IW(Wv,1)+1,1);
    T=wc./PP*PNS;

    GNZI=find(sum(WI,1)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
    T=T(:,GNZI);
    %% Rebel control
    RC=RC(GNZI);

    
    
    B=zeros(IW(Wv,2)-IW(Wv,1)+1,2);
    RCC=find(RC==1);
    NCC=find(RC==0);
    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_RebelControl_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Non-Rebel Control Govenorate Impacted \t Total Gov. Not impacted start \t Total NRC Gov. Not impacted start \t Prob. Impacted \t  Expected prob for NRCG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_RebelControl_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.Non_RebelControlGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForNRCG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,1)=1-chi2cdf(LR,1);
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(3,2,1);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['Rebel control (n=' num2str(length(RCC)) ')'],['Government control (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'RebelC_SA.png','-dpng','-r600');

    %% Conflict
    MC=mean(Ctv,2);
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));
CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Conflict_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low Conflict Govenorate Impacted \t Total Gov. Not impacted start \t Total LC Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LCG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Conflict_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.LowConflictGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLCG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,2)=1-chi2cdf(LR,1);

    subplot(3,2,2);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['High conflict (n=' num2str(length(RCC)) ')'],['Low conflict (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Conflict_SA.png','-dpng','-r600');

    %% Attacks
    MC=sum(tA,2);
    MC=MC(GNZI);

    RCC=find(MC>0);
    NCC=find(MC==0);

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Attack_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t No Attack Govenorate Impacted \t Total Gov. Not impacted start \t Total NA Gov. Not impacted start \t Prob. Impacted \t  Expected prob for NAG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Attack_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.NoAttackGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForNAG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,3)=1-chi2cdf(LR,1);

    subplot(3,2,3);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['Attacks (n=' num2str(length(RCC)) ')'],['No attacks (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Attack_SA.png','-dpng','-r600');

    %% Rainfall
    MC=mean(Rtv,2);
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Rain_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low Rain Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LRG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Rain_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.LowRainGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLRG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,4)=1-chi2cdf(LR,1);

    subplot(3,2,4);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['High rainall (n=' num2str(length(RCC)) ')'],['Low rainfall (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Rain_SA.png','-dpng','-r600');

    %% population density
    MC=P;
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    
    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_PD_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low PD Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LPDG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_PD_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.LowPDGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLPDG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,5)=1-chi2cdf(LR,1);

    subplot(3,2,5);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['High population density (n=' num2str(length(RCC)) ')'],['Low population density (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'PopD_SA.png','-dpng','-r600');


    %% health sites
    MC=H;
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Health_Wave-' num2str(Wv) 'CumulativeCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low HS Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LHSG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        f=find(T(ii,RCC)<=TES(Wv));
        B(ii,1)=length(f)./length(RCC);
        f=find(T(ii,NCC)<=TES(Wv));
        B(ii,2)=length(f)./length(NCC);
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Health_Wave-' num2str(Wv) 'CumulativeCases.dat']);
    O2=sum(SAT.LowHSGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLHSG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,6)=1-chi2cdf(LR,1);

    subplot(3,2,6);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' cumulative cases'],['per ' num2str(PNS)]})
    legend({['High no. health sites (n=' num2str(length(RCC)) ')'],['Low no. health sites (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    print(gcf,['SA-Wave_' num2str(Wv) '-Cumulative.png'],'-dpng','-r600');
end
fmf = fopen([pwd '\Tables\Pvalue-SA_Cumulative.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Wave \t Rebel_control \t Conflict \t Attacks \t Rain \t Pop_Density \t  Health sites \t \n'); % The header for the table
    fprintf(fmf, 'First \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(1,:)]); % The header for the table
    fprintf(fmf, 'Second \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(2,:)]); % The header for the table
    fprintf(fmf, 'Third \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(3,:)]); % The header for the table
    fprintf(fmf, 'Fourth \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t ',[PT(4,:)]); % The header for the table
    fprintf(fmf, 'All \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t ',[PT(5,:)]); % The header for the table
    fclose(fmf);
    
 readtable([pwd '\Tables\Pvalue-SA_Cumulative.dat'])
    


%% Weekly Cases

IW=[1 21; 22 74; 75 116; 117 153; 1 153];
PNS=10000;
TES=PNS.*[0.00007 0.0015 0.000015 0.00075 0.0015];
PT=zeros(5,6);
for Wv=1:5
    load('Yemen_Gov_Incidence.mat'); % Incidence data
    load('Yemen_National_Incidence.mat'); % Load the national incidence data for the comparison in the projection
    S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
    load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
    WI=IData; % Transpose the data set such that the number of areas is the row
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

    load('PopulationSize_Yemen.mat')
    load('Yemen_Gov_Incidence.mat')
    wc=(IData(IW(Wv,1):IW(Wv,2),:));
    PP=repmat(AP,IW(Wv,2)-IW(Wv,1)+1,1);
    T=wc./PP*PNS;

    GNZI=find(sum(WI,1)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
    T=T(:,GNZI);
    %% Rebel control
    RC=RC(GNZI);

    
    
    B=zeros(IW(Wv,2)-IW(Wv,1)+1,2);
    RCC=find(RC==1);
    NCC=find(RC==0);
    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_RebelControl_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Non-Rebel Control Govenorate Impacted \t Total Gov. Not impacted start \t Total NRC Gov. Not impacted start \t Prob. Impacted \t  Expected prob for NRCG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_RebelControl_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.Non_RebelControlGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForNRCG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,1)=1-chi2cdf(LR,1);
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(3,2,1);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['Rebel control (n=' num2str(length(RCC)) ')'],['Government control (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'RebelC_SA.png','-dpng','-r600');

    %% Conflict
    MC=mean(Ctv,2);
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));
CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Conflict_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low Conflict Govenorate Impacted \t Total Gov. Not impacted start \t Total LC Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LCG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Conflict_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.LowConflictGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLCG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,2)=1-chi2cdf(LR,1);

    subplot(3,2,2);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['High conflict (n=' num2str(length(RCC)) ')'],['Low conflict (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Conflict_SA.png','-dpng','-r600');

    %% Attacks
    MC=sum(tA,2);
    MC=MC(GNZI);

    RCC=find(MC>0);
    NCC=find(MC==0);

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Attack_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t No Attack Govenorate Impacted \t Total Gov. Not impacted start \t Total NA Gov. Not impacted start \t Prob. Impacted \t  Expected prob for NAG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Attack_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.NoAttackGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForNAG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,3)=1-chi2cdf(LR,1);

    subplot(3,2,3);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['Attack(s) (n=' num2str(length(RCC)) ')'],['No attacks (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Attack_SA.png','-dpng','-r600');

    %% Rainfall
    MC=mean(Rtv,2);
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Rain_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low Rain Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LRG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Rain_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.LowRainGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLRG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,4)=1-chi2cdf(LR,1);

    subplot(3,2,4);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['High rainfall (n=' num2str(length(RCC)) ')'],['Low rainfall (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'Rain_SA.png','-dpng','-r600');

    %% population density
    MC=P;
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    
    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_PD_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low PD Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LPDG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_PD_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.LowPDGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLPDG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,5)=1-chi2cdf(LR,1);

    subplot(3,2,5);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['High pop. density (n=' num2str(length(RCC)) ')'],['Low pop. density (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    %print(gcf,'PopD_SA.png','-dpng','-r600');


    %% health sites
    MC=H;
    MC=MC(GNZI);

    RCC=find(MC>=median(MC));
    NCC=find(MC<median(MC));

    CC=length(RCC)+length(NCC);
    CC2=length(NCC);
    Ttemp=T;
    fmf = fopen([pwd '\Tables\Survival_Health_Wave-' num2str(Wv) 'WeeklyCases.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Time of event \t Total Govnorate Impacted \t Low HS Govenorate Impacted \t Total Gov. Not impacted start \t Total LR Gov. Not impacted start \t Prob. Impacted \t  Expected prob for LHSG \t \n'); % The header for the table
    for ii=1:(IW(Wv,2)-IW(Wv,1)+1)
        f1=find(Ttemp(ii,RCC)>TES(Wv));
        f2=find(Ttemp(ii,NCC)>TES(Wv));
        Ttemp(:,RCC(f1))=0;
        Ttemp(:,NCC(f2))=0;
        
        if(~isempty(f1))
            fprintf(fmf, '%d \t %d \t %d \t %d \t %d \t %5.4f \t %5.4f \t\n',[IW(Wv,1)+(ii-1) length(f1)+length(f2) length(f2) CC CC2 (length(f1)+length(f2))./CC CC2.*(length(f1)+length(f2))./CC]); % The header for the table
            CC=CC-(length(f1)+length(f2));
            CC2=CC2-length(f2);
        end
        
        B(ii,1)=(CC-CC2)./length(RCC);
        B(ii,2)=CC2./length(NCC);
    end
    fclose(fmf);
    SAT=readtable([pwd '\Tables\Survival_Health_Wave-' num2str(Wv) 'WeeklyCases.dat']);
    O2=sum(SAT.LowHSGovenorateImpacted);
    O1=sum(SAT.TotalGovnorateImpacted)-O2;
    E2=sum(SAT.ExpectedProbForLHSG);
    E1=O1+O2-E2;
    LR=(O1-E1).^2./E1+(O2-E2)^2./E2;
    PT(Wv,6)=1-chi2cdf(LR,1);

    subplot(3,2,6);
    plot([IW(Wv,1):IW(Wv,2)],B(:,1),'r',[IW(Wv,1):IW(Wv,2)],B(:,2),'k','LineWidth',2); hold on
   % plot(75.*ones(101,1),[0:0.01:1],'b-.',22.*ones(101,1),[0:0.01:1],'b-.','LineWidth',1.5);
    set(gca,'Tickdir','out','LineWidth',2,'Fontsize',16);
    xlim([IW(Wv,1) IW(Wv,2)]);
    ylim([0 1]);
    box off
    xlabel('Week');
     ylabel({'Fraction of gov. with',[num2str(TES(Wv))  ' Weekly cases'],['per ' num2str(PNS)]})
    legend({['High no. health sites (n=' num2str(length(RCC)) ')'],['Low no. health sites (n=' num2str(length(NCC)) ')']},'Fontsize',16)
    legend boxoff
    print(gcf,['SA-Wave_' num2str(Wv) '-Weekly.png'],'-dpng','-r600');
end
fmf = fopen([pwd '\Tables\Pvalue-SA_Weekly.dat'],'w'); % Writing to the file name
    fprintf(fmf, 'Wave \t Rebel_control \t Conflict \t Attacks \t Rain \t Pop_Density \t  Health sites \t \n'); % The header for the table
    fprintf(fmf, 'First \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(1,:)]); % The header for the table
    fprintf(fmf, 'Second \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(2,:)]); % The header for the table
    fprintf(fmf, 'Third \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t \n',[PT(3,:)]); % The header for the table
    fprintf(fmf, 'All \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t  %5.4f \t ',[PT(4,:)]); % The header for the table
    fclose(fmf);
    
 readtable([pwd '\Tables\Pvalue-SA_Weekly.dat'])
    
%% Susceptiblility curve approximation

