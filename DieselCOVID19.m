function Dieselt=DieselCOVID19
%% Diesel prices for COVID-19
load('Diesel_Gov_Yemen_COVID-19.mat');
Diesel=Diesel'-min(Diesel(Diesel>0));
Dieselt=zeros(22,247);

% The dates for the prices being used
startData = [ ];
endData = [ ];
for yy=2016:2021
    if(yy==2016)
        for mm=10:12
            startData=[startData; datenum(yy,mm,1)];
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end
        end
    elseif(yy==2021)        
        for mm=1:6
            startData=[startData; datenum(yy,mm,1)];    
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end        
        end
    else
        for mm=1:12
            startData=[startData; datenum(yy,mm,1)];
            if(mm<12)
                endData=[endData; datenum(yy,mm+1,1)-1]; % We go to the first day of the next month and subtract one day
            else
                endData=[endData; datenum(yy+1,1,1)-1]; % We go to the first day of the next month and subtract one day
            end
        end
    end
end

EndFirstEpiWeek= datenum('10-09-2016');% Start of second epiweek is oct 10, 2016
for ii=1:247
    f=find(EndFirstEpiWeek+7*(ii-1)<=endData,1);% Need the first one to satisfy this condition as we increase over time and will elminate the other past months
    Dieselt(:,ii)=Diesel(:,f);
end
end