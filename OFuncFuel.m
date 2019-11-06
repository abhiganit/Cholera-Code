function F= OFuncFuel(x,DCt,DSt,Port,Strt,Ed,Ft)
beta=10.^(x(1:4));

CFt=zeros(size(Ft(:,1:end)));
SFt=zeros(size(Ft(:,1:end)));
Nt=zeros(size(Ft(:,1:end)));
for ii=1:length(Strt)
    CFt(:,ii)= DailyConflictFunction(Strt(ii):Ed(ii),DCt,0,0);
    SFt(:,ii)= DailyConflictFunction(Strt(ii):Ed(ii),DSt,0,0);
%     Nt(:,ii)= DailyConflictFunction(Strt(ii):Ed(ii),DSt,Func(3),XT(3));
end

PCt=repmat([0 0 0 0 0 0 0 0 0 0 0 0 0 exp(-([0:9 0:12])./(10.^(x(5))))],length(CFt(:,1)),1);

% for ii=1:length(Port)
%     f=find(Port(ii,:)>0);
%     if(length(f)==1)
%        PCt(ii,:)= Nt(f,:);
%     else
%        PCt(ii,:)= min(Nt(f,:)); 
%     end
% end

Y=beta(1).*ones(size(CFt))+beta(2).*CFt+beta(3).*SFt+beta(4).*PCt;

F=(Ft-Y).^2;
F=log10(sum(F(:)));
end
