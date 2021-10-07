clear;
close all;
rng default;
load('Fit-Vaccination-IncidenceperCapita-Conflict-Shellings-Diesel-Rain-CalibratedDAR.mat');
[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau,PopS,CI] = LoadYemenData;

peffect=zeros(4,1);
ptime=zeros(4,1);
medianeffect=zeros(4,3);
mediantime=zeros(4,3);
DateCOVID=cell(4,2);
DatePreCOVID=cell(4,2);
Cov=cell(4,1);
% Evaluate the number of paramters that are being used in the estimation 
[~,beta,tau,DB,DA,K,n,KP,KV,dV,r0,~,w]=RetParameterPS(par,XU,CF,maxtau);
% Conflict 


load('Area_Yemen.mat'); % 
load('PopulationSize_Yemen.mat'); % Populatino szie for 2016, 2017, 2018 and 2019 for the govneroates

PopS=AP(:,end);


S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weekly conflict
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Conflict_COVID-19_Timeline.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,258);

Ctv=log(1+Ctv./repmat(A,1,length(Ctv(1,:))).*PopS);

Ctv=Ctv(GNZI,:);

X=zeros(4,length(Ctv(:,1)),length(Ctv(1,(1+maxtau-tau(1)):(end-tau(1)))));
for ii=(maxtau+1):2*maxtau
    X(ii-maxtau,:,:)=beta(ii).*ImpactConflict(Ctv(:,(1+maxtau-tau(ii)):(end-tau(ii))),K(1),CF(1));
end

Xt=squeeze(X(1,:,:))+squeeze(X(2,:,:))+squeeze(X(3,:,:))+squeeze(X(4,:,:));
Yt=sum(Xt,2);
ff=find(Yt>0);
mtime=[];
for ii=1:length(ff)
    mtime=[mtime;[sum(Xt(ff(ii),[180:254])>0)]./sum(Xt(ff(ii),[105:179])>0)];
end

meffect=[];
for ii=1:length(ff)
    meffect=[meffect;[sum(Xt(ff(ii),[180:254]),2)]./sum(Xt(ff(ii),[105:179]),2)];
end

[peffect(1)]=signrank(meffect-1);
[ptime(1)]=signrank(mtime-1);
medianeffect(1,:)=[median(meffect) min(meffect) max(meffect)];
mediantime(1,:)=[median(mtime) min(mtime) max(mtime)];
DateCOVID(1,:)={'April 6, 2020','September 12, 2021'}; 
DatePreCOVID(1,:)={'October 29, 2018','April 5, 2020'}; 
Cov{1}={'Weekly conflict'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sheeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Yemen_Air_Shelling_COVID-19.mat');
Mt=GLevelConflict(YASt,S,258); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
Mt=log(1+Mt./repmat(A,1,length(Mt(1,:))).*PopS);

Mt=Mt(GNZI,:);
X2=zeros(4,length(Mt(:,1)),length(Mt(1,(1+maxtau-tau(1)):(end-tau(1)))));
for ii=(2.*maxtau+1):3*maxtau
    X2(ii-2.*maxtau,:,:)=beta(ii).*ImpactConflict(Mt(:,(1+maxtau-tau(ii)):(end-tau(ii))),K(2),CF(2));
end

X2t=squeeze(X2(1,:,:))+squeeze(X2(2,:,:))+squeeze(X2(3,:,:))+squeeze(X2(4,:,:));

Y2t=sum(X2t,2);
ff=find(Y2t>0);
mtime=[];

for ii=1:length(ff)
    mtime=[mtime;[sum(X2t(ff(ii),[180:254])>0)]./sum(X2t(ff(ii),[105:179])>0)];
end

meffect=[];
for ii=1:length(ff)
    meffect=[meffect;[sum(X2t(ff(ii),[180:254]),2)]./sum(X2t(ff(ii),[105:179]),2)];
end

[peffect(2)]=signrank(meffect-1);
[ptime(2)]=signrank(mtime-1);
medianeffect(2,:)=[median(meffect) min(meffect) max(meffect)];
mediantime(2,:)=[median(mtime) min(mtime) max(mtime)];
DateCOVID(2,:)={'April 6, 2020','September 12, 2021'}; 
DatePreCOVID(2,:)={'October 29, 2018','April 5, 2020'}; 
Cov{2}={'Shelling/attacks'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diesel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dieselt=DieselCOVID19;
Dieselt=Dieselt(GNZI,:);
X3=zeros(4,length(Dieselt(:,1)),length(Dieselt(1,(1+maxtau-tau(1)):(end-tau(1)))));
for ii=(3.*maxtau+1):4*maxtau
    X3(ii-3.*maxtau,:,:)=beta(ii).*max(Dieselt(:,(1+maxtau-tau(ii)):(end-tau(ii)))-KP(1),0);
end

X3t=squeeze(X3(1,:,:))+squeeze(X3(2,:,:))+squeeze(X3(3,:,:))+squeeze(X3(4,:,:));


Y3t=sum(X3t,2);
ff=find(Y3t>0);
mtime=[];

for ii=1:length(ff)
    mtime=[mtime;[sum(X3t(ff(ii),[180:243])>0)]./sum(X3t(ff(ii),[116:179])>0)];
end

meffect=[];
for ii=1:length(ff)
    meffect=[meffect;[sum(X3t(ff(ii),[180:243]),2)]./sum(X3t(ff(ii),[116:179]),2)];
end


[peffect(3)]=signrank(meffect-1);
[ptime(3)]=signrank(mtime-1);
medianeffect(3,:)=[median(meffect) min(meffect) max(meffect)];
mediantime(3,:)=[median(mtime) min(mtime) max(mtime)];
DateCOVID(3,:)={'April 6, 2020','June 27, 2021'}; 
DatePreCOVID(3,:)={'January 14, 2019','April 5, 2020'}; 
Cov{3}={'Diesel price'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rainfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Rainfall_COVID-19.mat');
Rtv=RainCOVID19';
Rtv=Rtv(GNZI,:);
X4=zeros(4,length(Rtv(:,1)),length(Rtv(1,(1+maxtau-tau(1)):(end-tau(1)))));
for ii=(5.*maxtau+1):6*maxtau
    X4(ii-5.*maxtau,:,:)=beta(ii).*ImpactRainfall(Rtv(:,(1+maxtau-tau(ii)):(end-tau(ii))),RF,r0);
end

X4t=squeeze(X4(1,:,:))+squeeze(X4(2,:,:))+squeeze(X4(3,:,:))+squeeze(X4(4,:,:));


Y4t=sum(X4t,2);
ff=find(Y4t>0);

mtime=[];

for ii=1:length(ff)
    mtime=[mtime;[sum(X4t(ff(ii),[180:217])>0)]./sum(X4t(ff(ii),[142:179])>0)];
end
meffect=[];
for ii=1:length(ff)
    meffect=[meffect;[sum(X4t(ff(ii),[180:217]),2)]./sum(X4t(ff(ii),[142:179]),2)];
end


[peffect(4)]=signrank(meffect-1);
[ptime(4)]=signrank(mtime-1);
medianeffect(4,:)=[median(meffect) min(meffect) max(meffect)];
mediantime(4,:)=[median(mtime) min(mtime) max(mtime)];
DateCOVID(4,:)={'April 6, 2020','December 27, 2020'}; 
DatePreCOVID(4,:)={'July 15, 2019','April 5, 2020'}; 
Cov{4}={'Rainfall'};

medianeffect=medianeffect-1;
mediantime=mediantime-1;

T=table(Cov,DatePreCOVID,DateCOVID,mediantime,ptime,medianeffect,peffect);
writetable(T,'COVID-19_Covariate_Changes.csv');
