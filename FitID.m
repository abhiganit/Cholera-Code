function F=FitID(x,indx,beta,tA,DB,DA,Ctv,K,n,tau,maxtau,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,V1,V2,KV,dV,Rtv,Temptv,r0,temp_0,WI,Pop,CI,DAR,w)

beta(indx)=10.^x;

   
[Yt,~]= LogisticModel(beta,tA,DB,DA,Ctv,K,n,tau,maxtau,WPIN,FPIN,Mt,Wheatt,Dieselt,KP,V1,V2,KV,dV,Rtv,Temptv,r0,temp_0,WI,Pop,CI,DAR,w);

F=(WI(:,(maxtau+1):end))-(Yt); % Compute the difference for the times and the locations that is tau weeks ahead
F=F(:);
end