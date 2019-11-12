function [lbga,ubga,lbps,ubps,IntC,part] = BoundsFitting(XU,parf)


flbps=[-32.*ones(1,length(XU)) zeros(1,length(XU)-7) -32.*ones(1,8) -32 -2 -32 -2 -32.*ones(1,7) -32 log10(0.9)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
fubps=[ 5.*ones(1,length(XU)) ones(1,length(XU)-7)  log10([ones(1,8) 20 1000 20 1000 120 120 120 120 120 1000 1 1000 exp(log(26/56)/(4*52))])]; % specify the upperbound for the parameters 

flb=[-32.*ones(1,length(XU)) ones(1,length(XU)-7)  -32.*ones(1,8) -32 -2 -32 -2 -32.*ones(1,7) -32 log10(0.9)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
fub=[ 5.*ones(1,length(XU)) 4.*ones(1,length(XU)-7)  log10([ones(1,8) 20 1000 20 1000 120 120 120 120 120 1000 1 1000 exp(log(26/56)/(4*52))])]; % specify the upperbound for the parameters 

IntCF=[1:((length(XU)-7))]+length(XU);
nob=length(XU);
Indx=zeros(size(fub));

Indx(1:length(XU))=XU;

Indx((1+length(XU)):(2.*length(XU)-7))=XU(8:end);

NTE=(length(XU)-7);
lenbeta=length(XU)+NTE;
if(sum(Indx(IntCF))>0)
    IntC=sum(XU)+[1:sum(Indx(IntCF))];
else
    IntC=[];
end
%% Attack asscoaited paramters
DA=zeros(2,1);
DB=zeros(2,1);
if(XU(8)==1)  
    Indx(lenbeta+[1 2])=1;
end

if(XU(14)==1)  
    
    Indx(lenbeta+[3 4])=1;
end

DAE=zeros(2,1);
DBE=zeros(2,1);

if(XU(10)==1)  % See if attacks being used at all
 
    Indx(lenbeta+[5 6])=1;

end

if(XU(15)==1)  % See if attacks being used at all

    Indx(lenbeta+[7 8])=1;
end



%% Conflict associated paramters
K=zeros(2,1);
n=zeros(2,1);
if(XU(9)==1) % See if conflict is being used at alls
        Indx(lenbeta+[9 10])=1;
end

if(XU(13)==1) % See if conflict is being used at alls
    
        Indx(lenbeta+[11 12])=1;
end


 %% Rainfall assocaited paramters 
rl=zeros(4,1);
if(XU(11)>=1) % See if rainfall is being used at all
    Indx(lenbeta+[13])=1;
end
if(XU(13)>=1) % See if rainfall is being used at all
    Indx(lenbeta+[14])=1; 
end
if(XU(14)>=1) % See if rainfall is being used at all
    Indx(lenbeta+[15])=1; 
end
if(XU(15)>=1) % See if rainfall is being used at all
    Indx(lenbeta+[16])=1;  
end

 %% Rainfall only
if(XU(12)>=1) % See if rainfall is being used at all
    Indx(lenbeta+[17])=1;
    
end

% Incidence per capita saturation
if(sum(XU([8:11 13:18]))>=1)
    Indx(lenbeta+[18])=1;
end

% Weight of WASH vs Food security
Indx(lenbeta+[19:21])=1;


lbga=flb(Indx==1);
ubga=fub(Indx==1);
lbps=flbps(Indx==1);
ubps=fubps(Indx==1);
part=parf(Indx==1);
end

