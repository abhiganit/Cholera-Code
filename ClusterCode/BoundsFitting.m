function [lbga,ubga,lbps,ubps,IntC,part] = BoundsFitting(XU,parf,CF,maxtau)


flbps=[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32];
fubps=[ 8.*ones(1,length(XU)) 0 log10(13.5) log10(152) log10(13.5) log10(152) log10(783) log10(400) log10(95.5) 5 log10(exp(log(26/56)/(4*52))) 0];

flb=[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 log10(0.5) -32];
fub=[ 8.*ones(1,length(XU)) 0 log10(13.5) log10(152.49) log10(13.5) log10(152.49) log10(783) log10(400) log10(95.5) 5 log10(exp(log(26/56)/(4*52))) 0];

Indx=zeros(size(fub));
Indx(1:length(XU))=XU;
IntC=[];
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

lbga=flb(Indx==1);
ubga=fub(Indx==1);
lbps=flbps(Indx==1);
ubps=fubps(Indx==1);
part=parf(Indx==1);
end

