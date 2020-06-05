function [par] = ExpandPar(part,XU,CF,maxtau,os)
% ExpandPar(part,XU,CF,maxtau,os) expands the parameter set based on the
% the truncated paramter set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% part - truncated parameter vector based on XU
% XU - vector of covairiates to include
% CF - conflict function being used
% maxtau - the maximum lag considered in the model
% os - optimization selection of ga or patternsearhc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% par - the full parameters set where the covariates not included have
% small value

if(os==1)
    par=[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32];
else
    par=[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32];
end

Indx=zeros(size(par));
Indx(1:length(XU))=XU;
lenbeta=length(XU);

%% Attack asscoaited paramters

if(sum(XU(1:maxtau))>=1)  % See if attacks being used at all
    Indx(lenbeta+1)=1;  %looking after the attack
     % Add estimated paramter
end



%% Conflict associated paramters
if(sum(XU((maxtau+1):2*maxtau))>=1) % See if conflict is being used at alls
    if(CF(1)~=0)
        Indx(lenbeta+2)=1;
    end
    if(CF(1)==2) % If the full hill function is being used
        Indx(lenbeta+3)=1; 
    end    
end

% Shellings
if(sum(XU((2.*maxtau+1):3*maxtau))>=1) % See if conflict is being used at alls
    if(CF(2)~=0)
        Indx(lenbeta+4)=1; 
    end
    if(CF(2)==2) % If the full hill function is being used
        Indx(lenbeta+5)=1; 
    end    
end


% Diesel price
if(sum(XU((3.*maxtau+1):4*maxtau))>=1)
    Indx(lenbeta+6)=1;
    
end

% Wheat price
if(sum(XU((4.*maxtau+1):5*maxtau))>=1)
    Indx(lenbeta+7)=1;
    
end

if(sum(XU((5.*maxtau+1):6*maxtau))>=1)
    Indx(lenbeta+[8])=1;
end
% Vaccination
Indx(lenbeta+[9])=1;
Indx(lenbeta+10)=1;

Indx(lenbeta+11)=1;

par(Indx==1)=part;

end

