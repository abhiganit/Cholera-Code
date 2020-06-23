close all;
clc;

[~,~,tA,~,~,~,RC,~,WPINm,FPINm,Dieselt,Wheatt,~,~,GNZI,~,~,~,~] = LoadYemenData;

S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen

load('Conflict_Yemen_Time_Location_Project_Forward.mat'); % Load the conflict in the area for the projection
Ctv=GLevelConflict(ProC,S,153);
load('Yemen_Air_Shelling.mat');
Mt=GLevelConflict(YASt,S,153); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize

% Need to increase the diesel price as we transformed it by subtracting the
% minimum
load('Diesel_Gov_Yemen.mat')
DDD=Diesel(:,GNZI)';
Dieselt=Dieselt+min(Diesel(Diesel>0));
load('Wheat_Gov_Yemen.mat')
WWW=Wheat(:,GNZI)';
Wheatt=Wheatt+min(Wheat(Wheat>0));

RC=RC(GNZI);
D=Dieselt(GNZI,:);
W=Wheatt(GNZI,:);
C=Ctv(GNZI,:);
S=Mt(GNZI,:);
Ws=WPINm(GNZI,:);
Fs=FPINm(GNZI,:);
TA=tA(GNZI,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Gov. level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


fprintf('========================================================================\n');
fprintf('Governorate level correlation \n');
fprintf('========================================================================\n');
% Conflict and Diesel

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(:,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(:,jj:end-(ii-jj));
    end
    DT=D(:,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Diesel and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Conflict and Wheat

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(:,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(:,jj:end-(ii-jj));
    end
    DT=W(:,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Wheat and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Diesel and Wheat

    T=DDD(:,1:end);
    DT=WWW(:,1:end);
    [rr,pp]=corr(T(:),DT(:),'Type','Kendall');
fprintf('Wheat and diesel (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

fprintf('========================================================================\n');
fprintf(' Rebel control Governorate level correlation \n');
fprintf('========================================================================\n');
% Conflict and Diesel

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(RC==1,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(RC==1,jj:end-(ii-jj));
    end
    DT=D(RC==1,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Diesel and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Conflict and Wheat

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(RC==1,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(RC==1,jj:end-(ii-jj));
    end
    DT=W(RC==1,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Wheat and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Diesel and Wheat


    T=DDD(RC==1,1:end);
    DT=WWW(RC==1,1:end);
    [rr,pp]=corr(T(:),DT(:),'Type','Kendall');
fprintf('Wheat and diesel (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

fprintf('========================================================================\n');
fprintf('Government controlled Governorate level correlation  \n');
fprintf('========================================================================\n');
% Conflict and Diesel

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(RC==0,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(RC==0,jj:end-(ii-jj));
    end
    DT=D(RC==0,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Diesel and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Conflict and Wheat

rt=zeros(70,1);
pt=zeros(70,1);
for ii=1:70
    T=C(RC==0,ii:end);
    for jj=(ii-1):-1:1
        T=T+C(RC==0,jj:end-(ii-jj));
    end
    DT=W(RC==0,ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end
fprintf('Wheat and conflict (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Diesel and Wheat


    T=DDD(RC==0,1:end);
    DT=WWW(RC==0,1:end);
    [rr,pp]=corr(T(:),DT(:),'Type','Kendall');
fprintf('Wheat and diesel (Governorate level): r=%3.2f , p=%3.2E \n',[rr pp]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Country. level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
fprintf('========================================================================\n');
fprintf('Country level correlation \n');
fprintf('========================================================================\n');
C=sum(C,1)';
S=sum(S,1)';
D=mean(D,1)';
W=mean(W,1)';
Ws=mean(Ws,1)';
Fs=mean(Fs,1)';
TA=sum(TA,1)';
% Diesel and conflict
for ii=1:70
    T=C(ii:end);
    for jj=(ii-1):-1:1
        T=T+C(jj:end-(ii-jj));
    end
    DT=D(ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end

fprintf('Diesel and conflict (Country level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Wheat and conflict
for ii=1:70
    T=C(ii:end);
    for jj=(ii-1):-1:1
        T=T+C(jj:end-(ii-jj));
    end
    DT=W(ii:end);
    [rt(ii),pt(ii)]=corr(T(:),DT(:),'Type','Kendall');
end
ff=find(pt<0.05);
if(~isempty(ff))    
   rr=max(rt(ff));
   gg=find(rr==rt(ff));
   pp=min(pt(ff(gg)));
else
   ff=find(pt==min(pt));
   pp=min(pt);
   rr=max(rt(ff));
end

fprintf('Wheat and conflict (Country level): r=%3.2f , p=%3.2E \n',[rr pp]);

% Wheat and diesle

    T=mean(DDD(:,1:end),1);
    DT=mean(WWW(:,1:end),1);
    [rr,pp]=corr(T(:),DT(:),'Type','Kendall');

fprintf('Wheat and diesel  (Country level): r=%3.2f , p=%3.2E \n',[rr pp]);