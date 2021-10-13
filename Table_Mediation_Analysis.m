clear;
load('Mediation_Analysis.mat');

CovariatePrice=CovariatePrice';
CovariateConflict=CovariateConflict';
ModelBase=ModelBase';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for correlation plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
betaz=cell(length(CovariateConflict),1);
betao=cell(length(CovariateConflict),1);
for ii=1:length(betaz)
    if(pValue_Corr(ii,1)<0.001)
        betaz{ii}={ [num2str(betaCorr(ii,1),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_Corr(ii,1)<0.01)
        betaz{ii}={ [num2str(betaCorr(ii,1),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_Corr(ii,1)<0.05)
        betaz{ii}={ [num2str(betaCorr(ii,1),'%5.4g') ,' \mbox{*}']};        
    else
        betaz{ii}={ [num2str(betaCorr(ii,1),'%5.4g')]};        
    end
        
    if(pValue_Corr(ii,2)<0.001)
        betao{ii}={ [num2str(betaCorr(ii,2),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_Corr(ii,2)<0.01)
        betao{ii}={ [num2str(betaCorr(ii,2),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_Corr(ii,2)<0.05)
        betao{ii}={ [num2str(betaCorr(ii,2),'%5.4g') ,' \mbox{*}']};        
    else
        betao{ii}={ [num2str(betaCorr(ii,2),'%5.4g')]};        
    end
end

 T=table(ModelBase,CovariatePrice,CovariateConflict,betaz,betao);
 
 writetable(T,'Mediation_Correlation.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betaone=cell(length(CovariateConflict),1);
betatwo=cell(length(CovariateConflict),1);
betathree=cell(length(CovariateConflict),1);
betafour=cell(length(CovariateConflict),1);

for ii=1:length(betaone)
    if(pValue_beta_Conflict(ii,1)<0.001)
        betaone{ii}={ [num2str(betaC(ii,1),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Conflict(ii,1)<0.01)
        betaone{ii}={ [num2str(betaC(ii,1),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Conflict(ii,1)<0.05)
        betaone{ii}={ [num2str(betaC(ii,1),'%5.4g') ,' \mbox{*}']};        
    else
        betaone{ii}={ [num2str(betaC(ii,1),'%5.4g')]};        
    end
    
    if(pValue_beta_Conflict(ii,2)<0.001)
        betatwo{ii}={ [num2str(betaC(ii,2),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Conflict(ii,2)<0.01)
        betatwo{ii}={ [num2str(betaC(ii,2),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Conflict(ii,2)<0.05)
        betatwo{ii}={ [num2str(betaC(ii,2),'%5.4g') ,' \mbox{*}']};        
    else
        betatwo{ii}={ [num2str(betaC(ii,2),'%5.4g')]};        
    end
    
    if(pValue_beta_Conflict(ii,3)<0.001)
        betathree{ii}={ [num2str(betaC(ii,3),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Conflict(ii,3)<0.01)
        betathree{ii}={ [num2str(betaC(ii,3),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Conflict(ii,3)<0.05)
        betathree{ii}={ [num2str(betaC(ii,3),'%5.4g') ,' \mbox{*}']};        
    else
        betathree{ii}={ [num2str(betaC(ii,3),'%5.4g')]};        
    end
    
    if(pValue_beta_Conflict(ii,4)<0.001)
        betafour{ii}={ [num2str(betaC(ii,4),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Conflict(ii,4)<0.01)
        betafour{ii}={ [num2str(betaC(ii,4),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Conflict(ii,4)<0.05)
        betafour{ii}={ [num2str(betaC(ii,4),'%5.4g') ,' \mbox{*}']};        
    else
        betafour{ii}={ [num2str(betaC(ii,4),'%5.4g')]};        
    end
end

T=table(ModelBase,CovariateConflict,betaone,betatwo,betathree,betafour);
 
 writetable(T,'Mediation_Conflict_Coefficients.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% construct symbols for the p-value for conflict auto-regression model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 
betaone=cell(length(CovariateConflict),1);
betatwo=cell(length(CovariateConflict),1);
betathree=cell(length(CovariateConflict),1);
betafour=cell(length(CovariateConflict),1);

for ii=1:length(betaone)
    if(pValue_beta_Mediation(ii,1)<0.001)
        betaone{ii}={ [num2str(betaC2(ii,1),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Mediation(ii,1)<0.01)
        betaone{ii}={ [num2str(betaC2(ii,1),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Mediation(ii,1)<0.05)
        betaone{ii}={ [num2str(betaC2(ii,1),'%5.4g') ,' \mbox{*}']};        
    else
        betaone{ii}={ [num2str(betaC2(ii,1),'%5.4g')]};        
    end
    
    if(pValue_beta_Mediation(ii,2)<0.001)
        betatwo{ii}={ [num2str(betaC2(ii,2),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Mediation(ii,2)<0.01)
        betatwo{ii}={ [num2str(betaC2(ii,2),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Mediation(ii,2)<0.05)
        betatwo{ii}={ [num2str(betaC2(ii,2),'%5.4g') ,' \mbox{*}']};        
    else
        betatwo{ii}={ [num2str(betaC2(ii,2),'%5.4g')]};        
    end
    
    if(pValue_beta_Mediation(ii,3)<0.001)
        betathree{ii}={ [num2str(betaC2(ii,3),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Mediation(ii,3)<0.01)
        betathree{ii}={ [num2str(betaC2(ii,3),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Mediation(ii,3)<0.05)
        betathree{ii}={ [num2str(betaC2(ii,3),'%5.4g') ,' \mbox{*}']};        
    else
        betathree{ii}={ [num2str(betaC2(ii,3),'%5.4g')]};        
    end
    
    if(pValue_beta_Mediation(ii,4)<0.001)
        betafour{ii}={ [num2str(betaC2(ii,4),'%5.4g') ,' \mbox{***}']};
    elseif(pValue_beta_Mediation(ii,4)<0.01)
        betafour{ii}={ [num2str(betaC2(ii,4),'%5.4g') ,' \mbox{**}']};
    elseif(pValue_beta_Mediation(ii,4)<0.05)
        betafour{ii}={ [num2str(betaC2(ii,4),'%5.4g') ,' \mbox{*}']};        
    else
        betafour{ii}={ [num2str(betaC2(ii,4),'%5.4g')]};        
    end
end

T=table(ModelBase,CovariatePrice,CovariateConflict,betaone,betatwo,betathree,betafour);
 
 writetable(T,'Mediation_Conflict_Comodity_Coefficients.csv','Delimiter',',');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 % Mediation effect
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
 
 EffectV=betaC-betaC2;
 pValue=zeros(8,4);
 Effectone=cell(8,1);
 Effecttwo=cell(8,1);
 Effectthree=cell(8,1);
 Effectfour=cell(8,1);
 for ii=1:length(Boot_betaC2)     
    MEv=Boot_betaC{ii}-Boot_betaC2{ii};
    ts=EffectV(ii,:)./std(MEv);
    pValue=1-tcdf(ts,length(XCv(:,1))-4);
    jj=1;
    if(pValue(jj)<0.001)
        Effectone{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}']};
    elseif(pValue(jj)<0.01)
        Effectone{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}']};
    elseif(pValue(jj)<0.05)
        Effectone{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}']};        
    else
        Effectone{ii}={ [num2str(EffectV(ii,jj),'%5.4g')]};        
    end
    
    jj=2;
    if(pValue(jj)<0.001)
        Effecttwo{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}']};
    elseif(pValue(jj)<0.01)
        Effecttwo{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}']};
    elseif(pValue(jj)<0.05)
        Effecttwo{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}']};        
    else
        Effecttwo{ii}={ [num2str(EffectV(ii,jj),'%5.4g')]};        
    end
    
    jj=3;
    if(pValue(jj)<0.001)
        Effectthree{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}']};
    elseif(pValue(jj)<0.01)
        Effectthree{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}']};
    elseif(pValue(jj)<0.05)
        Effectthree{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}']};        
    else
        Effectthree{ii}={ [num2str(EffectV(ii,jj),'%5.4g')]};        
    end
    
    jj=4;
    if(pValue(jj)<0.001)
        Effectfour{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{***}']};
    elseif(pValue(jj)<0.01)
        Effectfour{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{**}']};
    elseif(pValue(jj)<0.05)
        Effectfour{ii}={ [num2str(EffectV(ii,jj),'%5.4g') ,' \mbox{*}']};        
    else
        Effectfour{ii}={ [num2str(EffectV(ii,jj),'%5.4g')]};        
    end
 end
  T=table(ModelBase,CovariatePrice,CovariateConflict,Effectone,Effecttwo,Effectthree,Effectfour);
 
 writetable(T,'Mediation_Effect.csv','Delimiter',',');