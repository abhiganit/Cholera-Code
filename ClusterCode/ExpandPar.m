function [par] = ExpandPar(part,XU,os)
if(os==1)
    par=[-32.*ones(1,length(XU)) 10^(-32).*ones(1,length(XU)-7) 10^(-32).*ones(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,9) -32.*ones(1,3)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
else
    par=[-32.*ones(1,length(XU)) ones(1,length(XU)-7) zeros(1,7) -32.*ones(1,8) -32.*ones(1,4) -32.*ones(1,9) -32.*ones(1,3)]; 
end

nob=length(XU);
Indx=zeros(size(par));

Indx(1:length(XU))=XU;

Indx((1+length(XU)):(2.*length(XU)-7))=XU(8:end);


NTE=(length(XU)-7);
Indx(nob+NTE+[1 2])=XU([ 9 13]);
Indx(nob+NTE+3)=XU(12);
Indx(nob+NTE+[4 5 6 7])=XU([11 13 14 15]);
lenbeta=length(XU)+NTE+7;

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

% Threshold for Wheat
mln=zeros(3,1);
if(XU(16)>=1)
    Indx(lenbeta+[18])=1;
end
% Price threshold diesel
if(XU(17)>=1)
    Indx(lenbeta+[19])=1;
end
% Price threshold diesel
if(XU(18)>=1)
    Indx(lenbeta+[20])=1;
end

% Weight of WASH vs Food security
Indx(lenbeta+[21:24])=1;

par(Indx==1)=part;

end

