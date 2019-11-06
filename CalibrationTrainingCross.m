% Plots the error for the training and cross validation
close all;
atest=0;
NR=11;
CVEp=zeros(NR,1);
RSSp=zeros(NR,1);
for cvii=1:NR    
    PDS=0.4+0.05.*(cvii-1);
    load(['ForwardSelectionNoRainNoConflict-Vaccination-alpha=' num2str(atest*100) '-PercentData=' num2str(PDS*100) '.mat']);
    RSSp(cvii)=RSSv(end);
    CVEp(cvii)=CVE(end);
end

plot(0.4+[0:0.05:0.05.*(NR-1)],abs((RSSp+CVEp)-min(RSSp+CVEp))./2,'r','LineWidth',2); hold off;