function [lbga,ubga,IntC,part] = BoundsFitting(XU,parf,maxtau)
% BoundsFitting(XU,parf,maxtau) returns the bounds for the geneatic
% algorithm and pattern search based on the input vector XU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% XU - vector of covairiates to include
% parf - vector of all parameters in the model
% CF - conflict function being used
% maxtau - the maximum lag considered in the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%lbga - lowbound for genetic algiritm
%ubga - upperbound for genetic algiritm
%lbps - lowbound for pattern search
%ubps - uppoer bound for pattern search
%IntC - Integer contraints
%part - truncated paramters based on XU

flb=[-32.*ones(1,length(XU))  log10(0.25) -32 log10(0.5) -32 log10(0.5) -32 -32 -32 -32 -32 log10(0.5) -32 -6];
fub=[ 8.*ones(1,length(XU)) 0 log10(13.5) log10(152.49) log10(13.5) log10(152.49) log10(783) log10(400) log10(95.5) log10(40) 5 log10(exp(log(26/56)/(3*52))) 0 0];

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
   
        Indx(lenbeta+2)=1;
   
        Indx(lenbeta+3)=1; 
   
end

% Shellings
if(sum(XU((2.*maxtau+1):3*maxtau))>=1) % See if conflict is being used at alls
    
        Indx(lenbeta+4)=1; 
        Indx(lenbeta+5)=1; 
    
end


% Diesel price
if(sum(XU((3.*maxtau+1):4*maxtau))>=1)
    Indx(lenbeta+6)=1;
    
end

% Wheat price
if(sum(XU((4.*maxtau+1):5*maxtau))>=1)
    Indx(lenbeta+7)=1;
    
end

%Rainfall
Indx(lenbeta+[8])=1;


%Temprature 
Indx(lenbeta+[9])=1;

% Vaccination
Indx(lenbeta+[10])=1;
Indx(lenbeta+11)=1;


Indx(lenbeta+12)=1;
Indx(lenbeta+13)=1;

lbga=flb(Indx==1);
ubga=fub(Indx==1);
part=parf(Indx==1);
end

