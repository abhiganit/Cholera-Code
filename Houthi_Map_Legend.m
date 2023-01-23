plot(9.1+linspace(41.8,41.7741+0.8,1001),18.4.*ones(1001,1)-5,'k','LineWidth',2.5);
plot(9.1+linspace(41.8,41.7741+0.8,1001),18.95.*ones(1001,1)-5,'k','LineWidth',2.5);


plot(9.1+41.8.*ones(1001,1),linspace(18.4,18.95,1001)-5,'k','LineWidth',2.5);
plot(9.1+(41.7741+0.8).*ones(1001,1),linspace(18.4,18.95,1001)-5,'k','LineWidth',2.5);

plot(9.1+linspace(41.8,41.7741+0.8,1001),18.4+(18.95-18.4).*linspace(0,1,1001)-5,'k:','LineWidth',1.5);
plot(9.1+linspace(41.8,41.7741+0.8,1001),18.95-(18.95-18.4).*linspace(0,1,1001)-5,'k:','LineWidth',1.5);

text(9.1+41.7741+0.9,mean([18.4 18.95])-5,'Houthi Control','Fontsize',8);