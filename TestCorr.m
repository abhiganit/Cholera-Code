[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NW=length(WI(1,:));
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
close all;
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:64];
nd=WI(GNZI,(maxtau+1):end);
nd=length(nd(:));
CVE=10^6.*ones(length(INN),1);
AIC=10^6.*ones(length(INN),1);
MSS=10^6.*ones(length(INN),1);
k=zeros(length(INN),1);

[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDistrictData; % Load the data used to construct the figure
NW=length(WI(1,:));
load('Combo.mat');
for ii=1:length(INN)
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']))
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']);
        [lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par,CF,maxtau);
        MSS(ii)=RSSv;
        [k]=RetParameterPS(par,XU,CF,4);
        AIC(ii)=AICScore(k,nd,RSSv);
        CVE(ii)=OFuncDistrict(pars,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW));
    end
end
 AIC=AIC-min(AIC);
% FM=(MSS-mean(MSS))./std(abs(MSS));
FM=MSS;
% CVEv=(CVEv-mean(CVEv))./(std(CVEv));
ccc=[0.5 0.5 0.5;[221,28,119]./255;hex2rgb('#DE7A22');hex2rgb('#4C3F54');[153,52,4]./255; hex2rgb('#FAAF08');[5,112,176]./255;repmat([0 0 0],length(INN)-7,1)];
figure('units','normalized','outerposition',[0 0 1 1]);
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
for ii=1:length(INN)    
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']))
        scatter(FM(ii),CVE(ii),40,ccc(ii,:),'filled'); hold on    
        if(ii==1)
            text(FM(ii),CVE(ii),'-Incidence');
        else
            text(FM(ii),CVE(ii),[C(INC{INN(ii)}).N]);
        end
    end
end
box off;
set(gca,'LineWidth',2,'tickdir','out');
xlabel('MSE');
ylabel('CVE');

FM=AIC;
% CVEv=(CVEv-mean(CVEv))./(std(CVEv));
ccc=[0.5 0.5 0.5;[221,28,119]./255;hex2rgb('#DE7A22');hex2rgb('#4C3F54');[153,52,4]./255; hex2rgb('#FAAF08');[5,112,176]./255;repmat([0 0 0],length(INN)-7,1)];
figure('units','normalized','outerposition',[0 0 1 1]);
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
for ii=1:length(INN)    
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '.mat']))
        scatter(FM(ii),CVE(ii),40,ccc(ii,:),'filled'); hold on    
        if(ii==1)
            text(FM(ii),CVE(ii),'-Incidence');
        else
            text(FM(ii),CVE(ii),[C(INC{INN(ii)}).N]);
        end
    end
end
box off;
set(gca,'LineWidth',2,'tickdir','out');
xlabel('AIC');
ylabel('CVE');

MM=mean([MSS MSS],2);
TMM=sortrows([INN' MM],2);
figure('units','normalized','outerposition',[0 0 1 1]);
bar([1:length(MM)],TMM(:,2))
XTL={[C(INC{TMM(1,1)}).N]};
for ii=2:length(INN)
    if(isempty(INC{TMM(ii)}))
        XTL={XTL{1:end},'-Incidence'};
    else
        XTL={XTL{1:end},[C(INC{TMM(ii,1)}).N]};
    end
end
set(gca,'XTick',[1:length(MM)],'XTickLabel',XTL);
xtickangle(45);
