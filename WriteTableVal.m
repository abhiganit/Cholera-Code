%% AIC MSE CVE for the various models
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;
NWF=length(WI(1,:));
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par,CF,maxtau);

% Returns sum not average (will correct after adjusting the CVE)
MSE=OFuncProGAGLevel(pars,CF,WI(GNZI,1:NWF),tA(GNZI,1:NWF),Ctv(GNZI,1:NWF),XU,maxtau,WPIN(GNZI,1:NWF),FPIN(GNZI,1:NWF),Mt(GNZI,1:NWF),Wheatt(GNZI,1:NWF),Dieselt(GNZI,1:NWF),V1(GNZI,1:NWF),V2(GNZI,1:NWF),Rtv(GNZI,1:NWF),RF,PopS(GNZI,1:NWF),CI(GNZI,1:NWF));
MeanData_Fit=mean(WI(GNZI,1:NWF),2);


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenDataVal;
NW=length(WI(1,:));
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[lb,ub,lbps,ubps,IntC,pars] = BoundsFitting(XU,par,CF,maxtau);
CVE=OFuncProGAGLevel(pars,CF,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),XU,maxtau,WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW),Rtv(GNZI,1:NW),RF,PopS(GNZI,1:NW),CI(GNZI,1:NW));
MeanData_Val=mean(WI(GNZI,(NWF+1):NW),2);
CVE=CVE-MSE;

CVE=CVE./(NW-NWF);

MSE=MSE./(NWF-maxtau);

Norm_MSE=MSE./MeanData_Fit;
Norm_CVE=CVE./MeanData_Val;
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

Location={S(GNZI).ADM1_EN}';
 T=table(Location,MSE,MeanData_Fit,CVE,MeanData_Val,Norm_MSE,Norm_CVE);
 
 writetable(T,'TemporalValidation.csv','Delimiter',',');
