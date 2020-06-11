%% Runs the mediation analysis for diesel
clear;
close all;
clc;


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

load('Diesel_Gov_Yemen.mat')
Dieselt=Dieselt+min(Diesel(Diesel>0));

D=Dieselt(GNZI,:);
C=Ctv(GNZI,:);
S=Mt(GNZI,:);

D=D(:);
C=C(:);
S=S(:);


opts=optimset('MaxIter',10^4,'MaxFunEvals',10^4,'TolFun',10^(-9),'TolX',10^(-9),'UseParallel',true,'Display','off'); 
[bd,fval2]=fmincon(@(x)(mean((D-(x(1)+x(2)*C+x(3).*S)).^2)),([162 32 1]),[],[],[],[],[],[],[],opts);
bd=bd;
