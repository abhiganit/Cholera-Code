function [par] = ExpandPar(part,XU,os)
if(os==1)
    par=[-32.*ones(1,length(XU)) 10^(-32).*ones(1,length(XU)-7)  -32.*ones(1,8) -32 -4 -32 -4 -32.*ones(1,7)  -32 log10(0.9)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
else
    par=[-32.*ones(1,length(XU)) ones(1,length(XU)-7)  -32.*ones(1,8) -32 -4 -32 -4 -32.*ones(1,7) -32 log10(0.9)]; 
end

Indx=zeros(size(par));

Indx(1:length(XU))=XU;

Indx((1+length(XU)):(2.*length(XU)-7))=XU(8:end);


NTE=(length(XU)-7);
lenbeta=length(XU)+NTE;

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

par(Indx==1)=part;

end

