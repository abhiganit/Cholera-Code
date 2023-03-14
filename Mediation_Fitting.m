function F=Mediation_Fitting(x,D,C)
beta_j=10.^x;
Y=beta_j(1).*ones(size(C,2),size(C,3));
for jj=2:length(beta_j)
    Y=Y+beta_j(jj).*squeeze(C(jj-1,:,:));
end
Y=Y(:);
F=(D(:)-Y(:)).^2;

F=mean(F);
end