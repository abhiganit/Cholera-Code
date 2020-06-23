%% Produces table of AIC and Cross validatino Error


C=struct('N',{'-Targeted','-Conflict','-Shellings','-Diesel','-Wheat','-Rain'});
[WI,~,~,~,~,~,RC,~,~,~,~,~,~,~,GNZI,GV,maxtau,~,~] = LoadYemenData;
[GTF,GTCV] = SelectGov(WI,GNZI,GV,RC,0.8);
nd=WI(GNZI(GTF),(maxtau+1):end);
nd=length(nd(:));
AIC=zeros(4,4);
BIC=zeros(4,4);
CV=zeros(4,4);
indx=[2 3 5 6];
for yy=1:4
    for ss=yy:4
        load(['Fit-Vaccination-PercentData=80-IncidenceperCapita-Targeted-Diesel' C(unique([indx(ss) indx(yy)])).N '-Calibrate-DAR.mat']);
        [k]=RetParameterPS(par,XU,CF,4);
            CV(yy,ss)=CVE;
            CV(ss,yy)=CVE;
            AIC(yy,ss)=AICScore(k+1,nd,RSSv); % add one as DAR is estimated after
            AIC(ss,yy)=AICScore(k+1,nd,RSSv); % add one as DAR is estimated after
            BIC(yy,ss)=BICScore(k+1,nd,RSSv); % add one as DAR is estimated after
            BIC(ss,yy)=BICScore(k+1,nd,RSSv); % add one as DAR is estimated after
    end
end

clearvars -except AIC CV BIC indx
AIC=AIC-min(AIC(:));

f1=fopen('ModelCompTable-Next_CVEModel.txt','w');

C=struct('N',{'Targeted','Conflict','Shellings','Diesel','Wheat','Rain'});

fprintf(f1,'\\begin{sidewaystable}\\caption{Model comparison analysis using the Akaike information criterion (AIC) score and cross validation error (CVE)} \n \\begin{tabular}{c|cccccc} \n \\backslashbox{CVE}{$\\Delta$ AIC}');
for ii=1:4
    fprintf(f1,['& ' C(indx(ii)).N]);
end
fprintf(f1,'\\\\ \n \\hline ');

for yy=1:4
    fprintf(f1,C(indx(yy)).N);
    for ss=1:4
        if(yy==ss)
            if((AIC(yy,ss)==min(AIC(:)))||(CV(yy,ss)==min(CV(:)))) 
                fprintf(f1,['&\\cellcolor{blue!20} \\backslashbox{%3.2f}{%4.1f}'],[CV(yy,ss) AIC(yy,ss)]);
            else
                fprintf(f1,['& \\backslashbox{%3.2f}{%4.1f}'],[CV(yy,ss) AIC(yy,ss)]);
            end
        elseif(yy<ss)
            if((AIC(yy,ss)==min(AIC(:))))
                fprintf(f1,['& \\multicolumn{1}{r}{\\cellcolor{blue!20}{%4.1f}} '],[AIC(yy,ss)]);                 
            else
                fprintf(f1,['& \\multicolumn{1}{r}{%4.1f} '],[AIC(yy,ss)]);    
            end            
        else
            if((CV(yy,ss)==min(CV(:))))                
                fprintf(f1,['&  \\multicolumn{1}{l}{\\cellcolor{blue!20}{%3.2f}} '],[CV(yy,ss)]);
            else            
                fprintf(f1,['&  \\multicolumn{1}{l}{%3.2f} '],[CV(yy,ss)]);
            end
        end
    end
    if(yy<4)
        fprintf(f1,['\\\\ \n']);
    else
        fprintf(f1,['\n']);
    end
end
fprintf(f1,'\\end{tabular} \n \\end{sidewaystable}');
fclose(f1);