clear;
clc;

           
    
%% Names of the covariates
X=struct('N',{'WaSH','WaSH and Targeted Attacks','WaSH and Conflict','WaSH and Shelling','WaSH and Diesel','Food security','Food security and Conflict','Food security and Shellings','Food security and Diesel','Food security and Wheat'});
XU=[ones(1,10) 0 0];
flbps=[-32.*ones(1,length(XU)) -32 -32 -32 -6 -6 -6 -6 -32 -32 -32 log(0.9) -32 -32 -32 -32 -32];
fubps=[ ones(1,length(XU)) 0 0 0 3 5 3 5 4 4 5 log(exp(log(26/56)/(4*52))) 0 6 6 3 0];
temppar=[-0.0676715210877675,0.892606608496656,-16.7374095418067,-27.8356659108569,-17.4704811962739,-0.174499791393323,-28.4629616998684,-21.4266091791339,-21.8409078184253,-2.82669416070186,-30,-30,-0.0382099586294805,-3.98248176675222,-8.59668401843276,-6,-6,-6,-6,-27.1362273702645,-18.3276621316131,0.00912933824980688,-0.00368872675592762,-1.24987063985247,1.96251608269483,-1.28836484887011,-32,-32];

temppar=repmat(temppar,10,1).*(1+0.02.*(rand(10,length(temppar))-0.5));


[WI,Ctv,tA,Rtv,Mt,P,RC,H,WPIN,FPIN,Dieselt,Wheatt,V1,V2,GNZI,GV,maxtau] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,0.8);
NW=153; % Allow the model to fit the entire outbreak and cross validate among the govnerorates floor(153*PDS);

   %% Run the projection
   
      CF=[0 1 2 0 1 2 0 1 2;
          0 0 0 1 1 1 2 2 2];
      par=zeros(9,length(flbps));
      RSSv=zeros(9,1);
      CVE=zeros(9,1);
      RF=[-1; -1];
        for jj=1:9
           parfor ii=1:9
                [par(ii,:),RSSv(ii),CVE(ii)] =ProFittingGA(XU,0.8,CF(:,ii),RF,temppar);
           end
            for kk=1:9
               for ii=1:9                   
                   [~,~,~,~,~,part] = BoundsFitting(XU,par(ii,:),CF(:,kk),RF);
                   Rtemp=(OFuncProPS(part,CF(:,kk),WI(GNZI(GTF),1:NW),tA(GNZI(GTF),1:NW),Ctv(GNZI(GTF),1:NW),XU,maxtau,RC(GNZI(GTF)),WPIN(GNZI(GTF),1:NW),FPIN(GNZI(GTF),1:NW),Mt(GNZI(GTF),1:NW),Wheatt(GNZI(GTF),1:NW),Dieselt(GNZI(GTF),1:NW),V1(GNZI(GTF),1:NW),V2(GNZI(GTF),1:NW),Rtv(GNZI(GTF),1:NW),RF));
                   if(Rtemp<RSSv(kk))
                      RSSv(kk)=Rtemp; 
                      par(kk,:)=par(ii,:);
                      CVE(kk)=(OFuncProPS(part,CF(:,kk),WI(GNZI(GTCV),1:NW),tA(GNZI(GTCV),1:NW),Ctv(GNZI(GTCV),1:NW),XU,maxtau,RC(GNZI(GTCV)),WPIN(GNZI(GTCV),1:NW),FPIN(GNZI(GTCV),1:NW),Mt(GNZI(GTCV),1:NW),Wheatt(GNZI(GTCV),1:NW),Dieselt(GNZI(GTCV),1:NW),V1(GNZI(GTCV),1:NW),V2(GNZI(GTCV),1:NW),Rtv(GNZI(GTCV),1:NW),RF));
                   end
               end
            end
            temppar=unique(par,'rows');
            save(['Fit-Vaccination-PercentData=80_2.mat'],'par','RSSv','CVE','XU','X','CF','RF');
        end
