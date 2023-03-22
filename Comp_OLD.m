clear;


C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'});
 
Par_C=cell(2,32);
L=zeros(2,32);
NN=cell(1,32);
XUv=zeros(32,length(C));
load('Combo.mat');
for nn=1:32      
        indx=INC{nn};
        XUv(nn,indx)=1;
        load(['OLD-Fit-Vaccination-IncidenceperCapita' C(indx).N '-Rain-Calibrated_sigma.mat'])
        [~,~,~,part] = BoundsFitting(XU,par,4);
        L(1,nn)=-fval_c;
        load(['Fit-Vaccination-IncidenceperCapita' C(indx).N '.mat']);
        [~,~,~,part2] = BoundsFitting(XU,par,4);
        L(2,nn)=-fval_c;        
        Par_C{1,nn}=[part;part2];
        Par_C{2,nn}=1-part./part2;
        NN{nn}=[C(indx).N];
end

k=zeros(32,1);
for kk=1:32
k(kk)=length(Par_C{2,kk});
end

aics=aicbic(max(L,[],1),k);
daic=aics-min(aics);
w=exp(-daic./2)./sum(exp(-daic./2));
w_tot=zeros(1,length(C));
for nn=1:32
    w_tot=w_tot+w(nn).*XUv(nn,:);
end
clearvars -except Par_C L NN daic w_tot C


dL=L(1,:)-L(2,:);