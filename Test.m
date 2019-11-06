clear;
clc;
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenDistrictData;
PDS=0.90;
atest=0;
%% Forward selection
load(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
XU=XUr;
par=parr;

NW=length(WI(1,:));

TestE=(OFuncProPS(par,WI(1:21,1:NW),tA(1:21,1:NW),Ctv(1:21,1:NW),Rtv(1:21,1:NW),XU,maxtau,P(1:21,1:NW),RC(1:21),H(1:21,1:NW),WPIN(1:21,1:NW),FPIN(1:21,1:NW),Mt(1:21,1:NW),Wheatt(1:21,1:NW),Dieselt(1:21,1:NW),V1(1:21,1:NW),V2(1:21,1:NW)));

%% Load the data
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,PDS);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

TrainE=(OFuncProPS(par,WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),XU,maxtau,P(GNZI(GTF),1:NW),RC(GNZI(GTF)),H(GNZI(GTF),1:NW),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW)));
OVE=(OFuncProPS(par,WI(GNZI,1:NW),tA(GNZI,1:NW),Ctv(GNZI,1:NW),Rtv(GNZI,1:NW),XU,maxtau,P(GNZI,1:NW),RC(GNZI),H(GNZI,1:NW),WPIN(GNZI,1:NW),FPIN(GNZI,1:NW),Mt(GNZI,1:NW),Wheatt(GNZI,1:NW),Dieselt(GNZI,1:NW),V1(GNZI,1:NW),V2(GNZI,1:NW)));
CVE=(OFuncProPS(par,WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),XU,maxtau,P(GNZI(GTCV),1:NW),RC(GNZI(GTCV)),H(GNZI(GTCV),1:NW),WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW)));
bar([1 2 3 4],log10([TrainE CVE TestE OVE]))

