%% AIC MSE CVE for the various models
clear;
clc;
[WIF,Ctv,tA,Rtv,Temptv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZIF,GV,maxtau,PopS,CI] = LoadYemenData;
NW=length(WIF(1,:));

[WIVal,CtvVal,tAVal,RtvVal,TemptvVal,MtVal,PVal,RCVal,HVal,WPINVal,FPINVal,DieseltVal,WheattVal,V1Val,V2Val,GNZIVal,GVVal,maxtauVal,PopSVal,CIVal] = LoadYemenDataVal;
% The last two are governoerates that were used in the training of the
% model, we do not include them in the calcualtion of the cross validation
close all;
C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:32];
nd=WIF(GNZIF,(maxtau+1):end);
nd=length(nd(:));
MLLV=zeros(length(INN)+1,1);
AIC=zeros(length(INN)+1,1);
MLL=zeros(length(INN)+1,1);
Targeted=zeros(length(INN)+1,1);
Conflict=zeros(length(INN)+1,1);
Shellings=zeros(length(INN)+1,1);
Diesel=zeros(length(INN)+1,1);
Wheat=zeros(length(INN)+1,1);
k=zeros(length(INN)+1,1);
MN=zeros(length(INN)+1,1);

MN(end)=NaN;
MLL(end)=NaN;
MLLV(end)=NaN;

load('Combo.mat');
for ii=1:length(INN)
    if(isfile(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']))
        load(['Fit-Vaccination-IncidenceperCapita' C(INC{INN(ii)}).N '-CalibratedDAR.mat']);
        [lb,ub,IntC,pars] = BoundsFitting(XU,par,maxtau);
        
        
        tempD1=WIF(GNZIF,(maxtau+1):end);        
        
        tempD3=WIVal(GNZIVal,(154):end);  
        
        MLL(ii)=-RSSv./length(tempD1(:));
        
        [k(ii)]=RetParameterGA(par,XU,4);
        AIC(ii)=aicbic(-RSSv,k(ii));
        
        MLLV(ii)=-OFuncProGA_TempVal(pars,WIVal(GNZIVal,1:end),tAVal(GNZIVal,1:end),CtvVal(GNZIVal,1:end),XU,maxtau,WPINVal(GNZIVal,1:end),FPINVal(GNZIVal,1:end),MtVal(GNZIVal,1:end),WheattVal(GNZIVal,1:end),DieseltVal(GNZIVal,1:end),V1Val(GNZIVal,1:end),V2Val(GNZIVal,1:end),RtvVal(GNZIVal,1:end),TemptvVal(GNZIVal,1:end),PopSVal(GNZIVal,1:end),CIVal(GNZIVal,1:end));
        
        MLLV(ii)=MLLV(ii)./length(tempD3(:));
        
        if(ismember(1,INC{INN(ii)}))
            Targeted(ii)=1;
        else
            Targeted(ii)=0;
        end
        if(ismember(2,INC{INN(ii)}))
            Conflict(ii)=1;
        else
            Conflict(ii)=0;
        end
        if(ismember(3,INC{INN(ii)}))
            Shellings(ii)=1;
        else
            Shellings(ii)=0;
        end
        if(ismember(4,INC{INN(ii)}))
            Diesel(ii)=1;
        else
            Diesel(ii)=0;
        end
        if(ismember(5,INC{INN(ii)}))
            Wheat(ii)=1;
        else
            Wheat(ii)=0;
        end
        MN(ii)=ii;
    end
end
 AIC(1:end-1)=AIC(1:end-1)-min(AIC(1:end-1));
 AIC(end)=NaN;
 wAIC=exp(-AIC(1:end-1)./2)./sum(exp(-AIC(1:end-1)./2));
 wAIC=[wAIC;1];
 
Targeted(end)=sum(Targeted(1:end-1).*wAIC(1:end-1));
Conflict(end)=sum(Conflict(1:end-1).*wAIC(1:end-1));
Shellings(end)=sum(Shellings(1:end-1).*wAIC(1:end-1));
Diesel(end)=sum(Diesel(1:end-1).*wAIC(1:end-1));
Wheat(end)=sum(Wheat(1:end-1).*wAIC(1:end-1));

 T=table(MN,Targeted,Conflict,Shellings,Diesel,Wheat,MLL,MLLV,k,AIC,wAIC);
 
 writetable(T,'ModelFit.csv','Delimiter',',');
