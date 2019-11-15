function [lbga,ubga,lbps,ubps,IntC,part] = BoundsFitting(XU,parf,CF,RF)


flbps=[-32.*ones(1,length(XU)) -32 0 -1 0 -1 0 -1 0 -1 -32 -32 -32 log(0.9) -32 -32 -32 -32 -32];
fubps=[ 2.*ones(1,length(XU)) 0 log10(30) 1 log10(30) 1 log10(15) 1 log10(15) 1 log10(580) log10(300) 5 log(exp(log(26/56)/(4*52))) 0 6 6 log10(30) log10(30)];

flb=[-32.*ones(1,length(XU)) -32 0 -1 0 -1 0 -1 0 -1 -32 -32 -32 log(0.9) -32 -32 -32 -32 -32];
fub=[ 2.*ones(1,length(XU)) 0 log10(30) 1 log10(30) 1 log10(15) 1 log10(15) 1 log10(580) log10(300) 5 log(exp(log(26/56)/(4*52))) 0 6 6 log10(30) log10(30)];

IntC=[];
Indx=zeros(size(fub));
Indx(1:length(XU))=XU;
lenbeta=length(XU);
%% Attack asscoaited paramters
if(XU(2)==1)  % See if attacks being used at all
    Indx(lenbeta+1)=1;  %looking after the attack
end



%% Conflict associated paramters

if(XU(3)==1) % See if conflict is being used at alls
    if(CF(1)==2)
        Indx(lenbeta+2)=1;
    end
    if(CF(1)~=0) % If the full hill function is being used
        Indx(lenbeta+3)=1;
    end    
end

if(XU(7)==1) % See if conflict is being used at alls
    if(CF(2)==2)
        Indx(lenbeta+4)=1;
    end
    if(CF(2)~=0) % If the full hill function is being used
        Indx(lenbeta+5)=1; % Hill coefficient estimate
    end    
end
% Shellings
if(XU(4)==1)  % See if attacks being used at all
    if(CF(1)==2)
        Indx(lenbeta+6)=1; 
    end
    if(CF(1)~=0) % If the full hill function is being used
        Indx(lenbeta+7)=1; % Hill coefficient estimate
    end
end

if(XU(8)==1)   % See if attacks being used at all
    if(CF(2)==2)
        Indx(lenbeta+8)=1;
    end
    if(CF(2)~=0) % If the full hill function is being used
        Indx(lenbeta+9)=1; % Hill coefficient estimate
    end
end

KP=zeros(2,1);
% Diesel price
if(sum(XU([5 9])>=1))
   Indx(lenbeta+10)=1;
end

% Wheit price
if(XU(10)>=1)
   Indx(lenbeta+11)=1;
end

% Vaccination

Indx(lenbeta+[12])=1;
Indx(lenbeta+13)=1;

Indx(lenbeta+[14])=1;

if(RF(1)>=0)
    Indx(lenbeta+[15])=1;
end
if(RF(2)>=0)
    Indx(lenbeta+[16])=1;
end

lbga=flb(Indx==1);
ubga=fub(Indx==1);
lbps=flbps(Indx==1);
ubps=fubps(Indx==1);
part=parf(Indx==1);
end

