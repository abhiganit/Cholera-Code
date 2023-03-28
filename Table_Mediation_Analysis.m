clear;

load('Combo.mat');

% Models without Wheat but with Conflict or Shellings
ModelU={};
count=0;
ModelUW={};
for ii=1:length(INC)
    test=INC{ii};
    if(ismember(4,test) & (ismember(2,test)| ismember(3,test)))
        count=count+1;
        ModelUW{count}=[test];
        tt=[];
        for jj=1:length(test)
            if(test(jj)~=4)
                tt=[tt test(jj)];
            end
        end
        
        ModelU{count}=[tt];
    end
end

ModelNum=zeros(size(ModelU));
CovariateConflict=cell(size(ModelU));
CovariatePrice=cell(size(ModelU));

CName={'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat'};
CNameT={'Targeted/','Conflict/','Shellings/','Diesel/','Wheat/'};

load('Combo.mat');
for ii=1:length(ModelNum)
    CovariateConflict{ii}=[CNameT{[ModelU{ii}]} 'Incidence'];
    CovariatePrice{ii}='Diesel';
    for jj=1:length(INC)
        if(isequal(ModelU{ii},INC{jj}))
            ModelNum(ii)=jj;
        end
    end
end
    
CovariatePrice=CovariatePrice';
CovariateConflict=CovariateConflict';
ModelNum=ModelNum';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for correlation plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
beta_constant=cell(length(CovariateConflict),1);
beta_conflict=cell(length(CovariateConflict),1);
beta_shelling=cell(length(CovariateConflict),1);

for ii=1:length(beta_constant)
    load(['Mediation_Analysis' CName{ModelUW{ii}}  '.mat']);
    if(pValue_Corr(1)<0.001)
        beta_constant{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_Corr(1)<0.01)
        beta_constant{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_Corr(1)<0.05)
        beta_constant{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{*}'];        
    else
        beta_constant{ii}=[num2str(betaCorr(1),'%5.4g')];        
    end
    if(ismember(2,ModelUW{ii}))
        if(pValue_Corr(2)<0.001)
            beta_conflict{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{***}'];
        elseif(pValue_Corr(2)<0.01)
            beta_conflict{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{**}'];
        elseif(pValue_Corr(2)<0.05)
            beta_conflict{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{*}'];        
        else
            beta_conflict{ii}=[num2str(betaCorr(2),'%5.4g')];        
        end
    end
    if(ismember(3,ModelUW{ii}))
        if(pValue_Corr(end)<0.001)
            beta_shelling{ii}= [num2str(betaCorr(end),'%5.4g') ,' \mbox{***}'];
        elseif(pValue_Corr(end)<0.01)
            beta_shelling{ii}=[num2str(betaCorr(end),'%5.4g') ,' \mbox{**}'];
        elseif(pValue_Corr(end)<0.05)
            beta_shelling{ii}= [num2str(betaCorr(end),'%5.4g') ,' \mbox{*}'];        
        else
            beta_shelling{ii}=[num2str(betaCorr(end),'%5.4g')];        
        end
    end
end

 TC=table(ModelNum,CovariatePrice,CovariateConflict,beta_constant,beta_conflict,beta_shelling);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betacone=cell(length(CovariateConflict),1);
betactwo=cell(length(CovariateConflict),1);
betacthree=cell(length(CovariateConflict),1);
betacfour=cell(length(CovariateConflict),1);

betasone=cell(length(CovariateConflict),1);
betastwo=cell(length(CovariateConflict),1);
betasthree=cell(length(CovariateConflict),1);
betasfour=cell(length(CovariateConflict),1);

for ii=1:length(betacone)
    load(['Mediation_Analysis' CName{ModelUW{ii}}  '.mat'])
    if(pValue_beta_Conflict(1)<0.001)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(1)<0.01)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(1)<0.05)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{*}'];        
    else
        betacone{ii}=[num2str(betaC(1),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(2)<0.001)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(2)<0.01)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(2)<0.05)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{*}'];        
    else
        betactwo{ii}=[num2str(betaC(2),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(3)<0.001)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(3)<0.01)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(3)<0.05)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{*}'];        
    else
        betacthree{ii}=[num2str(betaC(3),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(4)<0.001)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(4)<0.01)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(4)<0.05)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{*}'];        
    else
        betacfour{ii}=[num2str(betaC(4),'%5.4g')];        
    end
    
    
    
    if(pValue_beta_Conflict(5)<0.001)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(5)<0.01)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(5)<0.05)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{*}'];        
    else
        betasone{ii}=[num2str(betaC(5),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(6)<0.001)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(6)<0.01)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(6)<0.05)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{*}'];        
    else
        betastwo{ii}=[num2str(betaC(6),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(7)<0.001)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(7)<0.01)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(7)<0.05)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{*}'];        
    else
        betasthree{ii}=[num2str(betaC(7),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(8)<0.001)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(8)<0.01)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(8)<0.05)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{*}'];        
    else
        betasfour{ii}=[num2str(betaC(8),'%5.4g')];        
    end
end

TM=table(ModelNum,CovariateConflict,betacone,betactwo,betacthree,betacfour,betasone,betastwo,betasthree,betasfour);
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betacone=cell(length(CovariateConflict),1);
betactwo=cell(length(CovariateConflict),1);
betacthree=cell(length(CovariateConflict),1);
betacfour=cell(length(CovariateConflict),1);


betasone=cell(length(CovariateConflict),1);
betastwo=cell(length(CovariateConflict),1);
betasthree=cell(length(CovariateConflict),1);
betasfour=cell(length(CovariateConflict),1);

for ii=1:length(betacone)
    load(['Mediation_Analysis' CName{ModelUW{ii}}  '.mat']);
    if(pValue_beta_Mediation(1)<0.001)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(1)<0.01)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(1)<0.05)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{*}'];        
    else
        betacone{ii}=[num2str(betaC2(1),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(2)<0.001)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(2)<0.01)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(2)<0.05)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{*}'];        
    else
        betactwo{ii}=[num2str(betaC2(2),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(3)<0.001)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(3)<0.01)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(3)<0.05)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{*}'];        
    else
        betacthree{ii}=[num2str(betaC2(3),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(4)<0.001)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(4)<0.01)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(4)<0.05)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{*}'];        
    else
        betacfour{ii}=[num2str(betaC2(4),'%5.4g')];        
    end
    
    
    
    if(pValue_beta_Mediation(5)<0.001)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(5)<0.01)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(5)<0.05)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{*}'];        
    else
        betasone{ii}=[num2str(betaC2(5),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(6)<0.001)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(6)<0.01)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(6)<0.05)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{*}'];        
    else
        betastwo{ii}=[num2str(betaC2(6),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(7)<0.001)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(7)<0.01)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(7)<0.05)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{*}'];        
    else
        betasthree{ii}=[num2str(betaC2(7),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(8)<0.001)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(8)<0.01)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(8)<0.05)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{*}'];        
    else
        betasfour{ii}=[num2str(betaC2(8),'%5.4g')];        
    end
end

TMC=table(ModelNum,CovariatePrice,CovariateConflict,betacone,betactwo,betacthree,betacfour,betasone,betastwo,betasthree,betasfour);
 
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 % Mediation effect
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 
 EffectV=betaC-betaC2;
 pValue=zeros(length(Boot_betaC2),8);
 Effectcone=cell(length(Boot_betaC2),1);
 Effectctwo=cell(length(Boot_betaC2),1);
 Effectcthree=cell(length(Boot_betaC2),1);
 Effectcfour=cell(length(Boot_betaC2),1);
 
 
 Effectsone=cell(length(Boot_betaC2),1);
 Effectstwo=cell(length(Boot_betaC2),1);
 Effectsthree=cell(length(Boot_betaC2),1);
 Effectsfour=cell(length(Boot_betaC2),1);
 
 for ii=1:length(Boot_betaC2)     
    MEv=Boot_betaC{ii}-Boot_betaC2{ii};
    ts=EffectV(ii,:)./std(MEv);
    pValue=1-tcdf(ts,length(XCv(:,1))-4);
    jj=1;
    if(pValue(jj)<0.001)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=2;
    if(pValue(jj)<0.001)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=3;
    if(pValue(jj)<0.001)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=4;
    if(pValue(jj)<0.001)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=5;
    if(pValue(jj)<0.001)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=6;
    if(pValue(jj)<0.001)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=7;
    if(pValue(jj)<0.001)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=8;
    if(pValue(jj)<0.001)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
 end
  TE=table(ModelNum,CovariatePrice,CovariateConflict,Effectcone,Effectctwo,Effectcthree,Effectcfour,Effectsone,Effectstwo,Effectsthree,Effectsfour);

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Wheat
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load('Mediation_Analysis_Alt_Wheat.mat');
ModelNum=zeros(size(ModelU));
CovariateConflict=cell(size(ModelU));
CovariatePrice=cell(size(ModelU));

CName={'Targeted-','Conflict-','Shellings-','Diesel-','Wheat-','Rain-'};

load('Combo.mat');
for ii=1:length(ModelNum)
    CovariateConflict{ii}=[CName{[ModelU{ii}]} 'Incidence'];
    CovariatePrice{ii}='Wheat';
    for jj=1:length(INC)
        if(isequal(ModelU{ii},INC{jj}))
            ModelNum(ii)=jj;
        end
    end
end
    
CovariatePrice=CovariatePrice';
CovariateConflict=CovariateConflict';
ModelNum=ModelNum';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for correlation plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
betacz=cell(length(CovariateConflict),1);
betaco=cell(length(CovariateConflict),1);

betasz=cell(length(CovariateConflict),1);
betaso=cell(length(CovariateConflict),1);

for ii=1:length(betacz)
    if(pValue_Corr(1)<0.001)
        betacz{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_Corr(1)<0.01)
        betacz{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_Corr(1)<0.05)
        betacz{ii}=[num2str(betaCorr(1),'%5.4g') ,' \mbox{*}'];        
    else
        betacz{ii}=[num2str(betaCorr(1),'%5.4g')];        
    end
        
    if(pValue_Corr(2)<0.001)
        betaco{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_Corr(2)<0.01)
        betaco{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_Corr(2)<0.05)
        betaco{ii}=[num2str(betaCorr(2),'%5.4g') ,' \mbox{*}'];        
    else
        betaco{ii}=[num2str(betaCorr(2),'%5.4g')];        
    end
    
    if(pValue_Corr(3)<0.001)
        betasz{ii}=[num2str(betaCorr(3),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_Corr(3)<0.01)
        betasz{ii}=[num2str(betaCorr(3),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_Corr(3)<0.05)
        betasz{ii}=[num2str(betaCorr(3),'%5.4g') ,' \mbox{*}'];        
    else
        betasz{ii}=[num2str(betaCorr(3),'%5.4g')];        
    end
        
    if(pValue_Corr(4)<0.001)
        betaso{ii}=[num2str(betaCorr(4),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_Corr(4)<0.01)
        betaso{ii}=[num2str(betaCorr(4),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_Corr(4)<0.05)
        betaso{ii}=[num2str(betaCorr(4),'%5.4g') ,' \mbox{*}'];        
    else
        betaso{ii}=[num2str(betaCorr(4),'%5.4g')];        
    end
end

 TCt=table(ModelNum,CovariatePrice,CovariateConflict,betacz,betaco,betasz,betaso);
 
 TC=[TC;TCt];
 writetable(TC,'Mediation_Correlation_Alt.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betacone=cell(length(CovariateConflict),1);
betactwo=cell(length(CovariateConflict),1);
betacthree=cell(length(CovariateConflict),1);
betacfour=cell(length(CovariateConflict),1);

betasone=cell(length(CovariateConflict),1);
betastwo=cell(length(CovariateConflict),1);
betasthree=cell(length(CovariateConflict),1);
betasfour=cell(length(CovariateConflict),1);

for ii=1:length(betacone)
    if(pValue_beta_Conflict(1)<0.001)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(1)<0.01)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(1)<0.05)
        betacone{ii}=[num2str(betaC(1),'%5.4g') ,' \mbox{*}'];        
    else
        betacone{ii}=[num2str(betaC(1),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(2)<0.001)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(2)<0.01)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(2)<0.05)
        betactwo{ii}=[num2str(betaC(2),'%5.4g') ,' \mbox{*}'];        
    else
        betactwo{ii}=[num2str(betaC(2),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(3)<0.001)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(3)<0.01)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(3)<0.05)
        betacthree{ii}=[num2str(betaC(3),'%5.4g') ,' \mbox{*}'];        
    else
        betacthree{ii}=[num2str(betaC(3),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(4)<0.001)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(4)<0.01)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(4)<0.05)
        betacfour{ii}=[num2str(betaC(4),'%5.4g') ,' \mbox{*}'];        
    else
        betacfour{ii}=[num2str(betaC(4),'%5.4g')];        
    end
    
    
    
    if(pValue_beta_Conflict(5)<0.001)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(5)<0.01)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(5)<0.05)
        betasone{ii}=[num2str(betaC(5),'%5.4g') ,' \mbox{*}'];        
    else
        betasone{ii}=[num2str(betaC(5),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(6)<0.001)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(6)<0.01)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(6)<0.05)
        betastwo{ii}=[num2str(betaC(6),'%5.4g') ,' \mbox{*}'];        
    else
        betastwo{ii}=[num2str(betaC(6),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(7)<0.001)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(7)<0.01)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(7)<0.05)
        betasthree{ii}=[num2str(betaC(7),'%5.4g') ,' \mbox{*}'];        
    else
        betasthree{ii}=[num2str(betaC(7),'%5.4g')];        
    end
    
    if(pValue_beta_Conflict(8)<0.001)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Conflict(8)<0.01)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Conflict(8)<0.05)
        betasfour{ii}=[num2str(betaC(8),'%5.4g') ,' \mbox{*}'];        
    else
        betasfour{ii}=[num2str(betaC(8),'%5.4g')];        
    end
end

TMt=table(ModelNum,CovariateConflict,betacone,betactwo,betacthree,betacfour,betasone,betastwo,betasthree,betasfour);
 TM=[TM;TMt];
 
 [C,IA,IC] = unique(TM.ModelNum);
 TM=TM(IA,:);
 writetable(TM,'Mediation_Conflict_Coefficients_Alt.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betacone=cell(length(CovariateConflict),1);
betactwo=cell(length(CovariateConflict),1);
betacthree=cell(length(CovariateConflict),1);
betacfour=cell(length(CovariateConflict),1);


betasone=cell(length(CovariateConflict),1);
betastwo=cell(length(CovariateConflict),1);
betasthree=cell(length(CovariateConflict),1);
betasfour=cell(length(CovariateConflict),1);

for ii=1:length(betacone)
    if(pValue_beta_Mediation(1)<0.001)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(1)<0.01)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(1)<0.05)
        betacone{ii}=[num2str(betaC2(1),'%5.4g') ,' \mbox{*}'];        
    else
        betacone{ii}=[num2str(betaC2(1),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(2)<0.001)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(2)<0.01)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(2)<0.05)
        betactwo{ii}=[num2str(betaC2(2),'%5.4g') ,' \mbox{*}'];        
    else
        betactwo{ii}=[num2str(betaC2(2),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(3)<0.001)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(3)<0.01)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(3)<0.05)
        betacthree{ii}=[num2str(betaC2(3),'%5.4g') ,' \mbox{*}'];        
    else
        betacthree{ii}=[num2str(betaC2(3),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(4)<0.001)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(4)<0.01)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(4)<0.05)
        betacfour{ii}=[num2str(betaC2(4),'%5.4g') ,' \mbox{*}'];        
    else
        betacfour{ii}=[num2str(betaC2(4),'%5.4g')];        
    end
    
    
    
    if(pValue_beta_Mediation(5)<0.001)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(5)<0.01)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(5)<0.05)
        betasone{ii}=[num2str(betaC2(5),'%5.4g') ,' \mbox{*}'];        
    else
        betasone{ii}=[num2str(betaC2(5),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(6)<0.001)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(6)<0.01)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(6)<0.05)
        betastwo{ii}=[num2str(betaC2(6),'%5.4g') ,' \mbox{*}'];        
    else
        betastwo{ii}=[num2str(betaC2(6),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(7)<0.001)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(7)<0.01)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(7)<0.05)
        betasthree{ii}=[num2str(betaC2(7),'%5.4g') ,' \mbox{*}'];        
    else
        betasthree{ii}=[num2str(betaC2(7),'%5.4g')];        
    end
    
    if(pValue_beta_Mediation(8)<0.001)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{***}'];
    elseif(pValue_beta_Mediation(8)<0.01)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{**}'];
    elseif(pValue_beta_Mediation(8)<0.05)
        betasfour{ii}=[num2str(betaC2(8),'%5.4g') ,' \mbox{*}'];        
    else
        betasfour{ii}=[num2str(betaC2(8),'%5.4g')];        
    end
end

TMCt=table(ModelNum,CovariatePrice,CovariateConflict,betacone,betactwo,betacthree,betacfour,betasone,betastwo,betasthree,betasfour);
 TMC=[TMC;TMCt];
 writetable(TMC,'Mediation_Conflict_Comodity_Coefficients_Alt.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 % Mediation effect
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 
 EffectV=betaC-betaC2;
 pValue=zeros(length(Boot_betaC2),8);
 Effectcone=cell(length(Boot_betaC2),1);
 Effectctwo=cell(length(Boot_betaC2),1);
 Effectcthree=cell(length(Boot_betaC2),1);
 Effectcfour=cell(length(Boot_betaC2),1);
 
 
 Effectsone=cell(length(Boot_betaC2),1);
 Effectstwo=cell(length(Boot_betaC2),1);
 Effectsthree=cell(length(Boot_betaC2),1);
 Effectsfour=cell(length(Boot_betaC2),1);
 
 for ii=1:length(Boot_betaC2)     
    MEv=Boot_betaC{ii}-Boot_betaC2{ii};
    ts=EffectV(ii,:)./std(MEv);
    pValue=1-tcdf(ts,length(XCv(:,1))-4);
    jj=1;
    if(pValue(jj)<0.001)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcone{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=2;
    if(pValue(jj)<0.001)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectctwo{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=3;
    if(pValue(jj)<0.001)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcthree{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=4;
    if(pValue(jj)<0.001)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectcfour{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=5;
    if(pValue(jj)<0.001)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsone{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=6;
    if(pValue(jj)<0.001)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectstwo{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=7;
    if(pValue(jj)<0.001)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsthree{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
    
    jj=8;
    if(pValue(jj)<0.001)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}'];
    elseif(pValue(jj)<0.01)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}'];
    elseif(pValue(jj)<0.05)
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}'];        
    else
        Effectsfour{ii}=[num2str(EffectV(ii,jj),'%5.4g')];        
    end
 end
  TEt=table(ModelNum,CovariatePrice,CovariateConflict,Effectcone,Effectctwo,Effectcthree,Effectcfour,Effectsone,Effectstwo,Effectsthree,Effectsfour);
 TE=[TE;TEt];
 writetable(TE,'Mediation_Effect_Alt.csv','Delimiter',',');