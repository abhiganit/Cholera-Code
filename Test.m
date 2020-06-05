C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
INN=[1:64];
load('Combo.mat');
ll=zeros(1,64);
for ii=1:64
   ll(ii)=length(INC{ii}); 
end
 fprintf(['=============================================================\n']);
 fprintf(['=============================================================\n']);
for ii=length(INN):-1:2    
    load(['Fit-Vaccination-IncidenceperCapita' C(INC{ii}).N '.mat'],'RSSv');
    MSE=RSSv;
    temp=INC{ii};
    f=combnk(temp,ll(ii)-1);
    for jj=1:length(f)
        load(['Fit-Vaccination-IncidenceperCapita' C(f(jj,:)).N '.mat'],'RSSv');
        if(MSE>RSSv)
            fprintf([C(INC{ii}).N ' --> ' C(f(jj,:)).N  ' %7.6f \n'],[(MSE-RSSv)]);
        end
    end
end