clear;
clc;

[WIa,Ctva,tAa,Rtva,Mta,Pa,RCa,Ha,WPINma,IDPta,GNZIa,maxtaua] = LoadYemenDistrictData;
HCity=[29 31 71];
%% Shellings from Janaury 1, 2015 to October 2, 2016


fprintf('===================================================== \n');
fprintf('Prior to epidemic \n');
fprintf('===================================================== \n');
load('PriorShelling_Yemen.mat');

% Adjust the timing index for the Prior shelling as they are currently
% negative
ShellP(:,1)=ShellP(:,1)-min(ShellP(:,1))+1;


S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen
% Sana'a City
fS=zeros(length(S),1);
for ii=1:length(fS)
  fS(ii)=strcmp({'Amanat Al Asimah'},S(ii).ADM1_EN); 
end
M=GLevelConflict(ShellP,S(fS==1),max(ShellP(:,1)));
fprintf(' Number of attacks from Janaury 1, 2015 to October 2, 2016 in Sanaa City: %d \n', sum(M(:)));
fprintf(' \n \n');
% Hodeidah City
M=GLevelConflict(ShellP,S(HCity),max(ShellP(:,1)));
fprintf(' Number of attacks from Janaury 1, 2015 to October 2, 2016 in Hodeidah City: %d \n', sum(M(:)));
fprintf(' \n \n');
% Aden
fA=zeros(length(S),1);
for ii=1:length(fA)
    fA(ii)=strcmp({'Aden'},S(ii).ADM1_EN);
end
M=GLevelConflict(ShellP,S(fA==1),max(ShellP(:,1)));
fprintf(' Number of attacks from Janaury 1, 2015 to October 2, 2016 in Aden: %d \n', sum(M(:)));
fprintf(' \n \n');
%% Shellings After October 3, 2016
load('Yemen_Air_Shelling.mat');

%% Shellings October 3, 2016 to April 9, 2017
fprintf('===================================================== \n');
fprintf('First wave \n');
fprintf('===================================================== \n');
startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-09-2017');% End date

% Sana'a City
fS=zeros(length(S),1);
for ii=1:length(fS)
  fS(ii)=strcmp({'Amanat Al Asimah'},S(ii).ADM1_EN); 
end
M=GLevelConflict(YASt,S(fS==1),ceil((1+endDateofSim-startDateofSim)./7)); 

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Sanaa City: %d \n'], sum(M(:)));

endDateofSim = datenum('11-20-2016');% End date

MP=GLevelConflict(YASt,S(fS==1),ceil((1+endDateofSim-startDateofSim)./7)); 
fprintf([' Percentage of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Sanaa City: %3.1f %% \n'], 100.*sum(MP(:))./sum(M(:)));

% Hodeadah City
fprintf(' \n \n');
startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-09-2017');% End date
fH=zeros(length(S),1);
for ii=1:length(fH)
  fH(ii)=strcmp({'Al Hudaydah'},S(ii).ADM1_EN); 
end

MH=GLevelConflict(YASt,S(fH==1),ceil((1+endDateofSim-startDateofSim)./7));

M=GLevelConflict(YASt,S(HCity),ceil((1+endDateofSim-startDateofSim)./7));

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Hodeadah City: %d \n'], sum(M(:)));
fprintf([' Percentage of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Hodeadah City relative to Governorate: %3.1f %% \n'], 100.*sum(M(:))./sum(MH(:)));
fprintf(' \n \n');
% Aden
fA=zeros(length(S),1);
for ii=1:length(fA)
    fA(ii)=strcmp({'Aden'},S(ii).ADM1_EN);
end

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-09-2017');% End date
M=GLevelConflict(YASt,S(fA==1),ceil((1+endDateofSim-startDateofSim)./7)); 

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden City: %d \n'], sum(M(:)));

endDateofSim = datenum('11-20-2016');% End date

MP=GLevelConflict(YASt,S(fA==1),ceil((1+endDateofSim-startDateofSim)./7)); 
fprintf([' Percentage of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden: %3.1f %% \n'], 100.*sum(MP(:))./sum(M(:)));
fprintf(' \n \n');
fprintf('===================================================== \n');
fprintf('Second wave \n');
fprintf('===================================================== \n');

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-23-2017');% End date

% Sana'a City
fS=zeros(length(S),1);
for ii=1:length(fS)
  fS(ii)=strcmp({'Amanat Al Asimah'},S(ii).ADM1_EN); 
end
M=GLevelConflict(YASt,S(fS==1),ceil((1+endDateofSim-startDateofSim)./7)); 

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Sanaa City: %d \n'], sum(M(:)));
temp=S(fS==1);
for ii=1:sum(fS)
    fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %d \n'], sum(M(ii,:)));
end

% May 1
startDateofSim = datenum('5-01-2017');% Start date
endDateofSim = datenum('4-29-2018');% End date
temp=Mta(23,1:ceil((1+endDateofSim-startDateofSim)./7));
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack from  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Sanaa City: %3.1f \n'],length(temp)./length(f));
fprintf(['Avg number of attacks from  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Sanaa City: %3.1f \n'],mean(temp));
temp=S(fS==1);
for ii=1:sum(fS)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii+3,1:ceil((1+endDateofSim-startDateofSim)./7))));
end

fprintf(' \n \n');
% Hodeadah City

startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-23-2017');% End date

M=GLevelConflict(YASt,S(HCity),ceil((1+endDateofSim-startDateofSim)./7)); 

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Hodeadah City: %d \n'], sum(M(:)));
temp=S(HCity);
for ii=1:length(HCity)
    fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %d \n'], sum(M(ii,:)));
end

% May 1
startDateofSim = datenum('5-01-2017');% Start date
endDateofSim = datenum('4-29-2018');% End date
temp=Mta(22,1:ceil((1+endDateofSim-startDateofSim)./7));
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Hodeadah City: %3.1f \n'],length(temp)./length(f));
fprintf(['Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Hodeadah City: %d \n'],sum(temp));
fprintf(['Avg number of attacks from  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Hodeadah City: %3.1f \n'],mean(temp));
temp=S(HCity);
for ii=1:length(HCity)
    fprintf(['Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %d \n'], sum(Mta(ii,1:ceil((1+endDateofSim-startDateofSim)./7))));
end
temp=S(HCity);
for ii=1:length(HCity)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii,1:ceil((1+endDateofSim-startDateofSim)./7))));
end
fprintf(' \n \n');

% Aden


startDateofSim = datenum('10-03-2016');% Start date
endDateofSim = datenum('4-23-2017');% End date

fA=zeros(length(S),1);
for ii=1:length(fA)
  fA(ii)=strcmp({'Aden'},S(ii).ADM1_EN); 
end
M=GLevelConflict(YASt,S(fA==1),ceil((1+endDateofSim-startDateofSim)./7)); 

fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden: %d \n'], sum(M(:)));
temp=S(fA==1);
for ii=1:sum(fA)
    fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %d \n'], sum(M(ii,:)));
end

% May 1


startDateofSim = datenum('5-01-2017');% Start date
endDateofSim = datenum('4-29-2018');% End date
fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden: %d \n'], sum(Mta(24,1:ceil((1+endDateofSim-startDateofSim)./7))));
temp=S(fA==1);
for ii=1:sum(fA)
    fprintf([' Number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %d \n'], sum(Mta(ii+3+sum(fS),1:ceil((1+endDateofSim-startDateofSim)./7))));
end

temp=Mta(24,1:ceil((1+endDateofSim-startDateofSim)./7));
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack from  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden: %3.1f \n'],length(temp)./length(f));
fprintf(['Avg number of attacks from  ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in Aden: %3.1f \n'],mean(temp));
temp=S(fA==1);
for ii=1:sum(fA)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii+3+sum(fS),1:ceil((1+endDateofSim-startDateofSim)./7))));
end

fprintf(' \n \n');


fprintf('===================================================== \n');
fprintf('Third and fourth wave \n');
fprintf('===================================================== \n');

startDateofSimA = datenum('5-01-2017');% Start date
endDateofSimA = datenum('4-29-2018');% End date

StartIndx=ceil((1+endDateofSimA-startDateofSimA)./7)+1;

startDateofSim = datenum('4-30-2018');% start date
endDateofSim = datenum('9-8-2019');% End date
EndIndx=ceil((1+endDateofSim-startDateofSimA)./7);
temp=Mta(23,(StartIndx):EndIndx);
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Sanaa City: %3.1f \n'],length(temp)./length(f));
fprintf(['Avg. number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Sanaa City: %3.1f \n'],mean(temp));
temp=S(fS==1);
for ii=1:sum(fS)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii+3,(StartIndx):EndIndx)));
end

fprintf(' \n \n');

% Hodeadah City

temp=Mta(22,(StartIndx):EndIndx);
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Hodeadah City: %3.1f \n'],length(temp)./length(f));
fprintf(['Avg. number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Hodeadah City: %3.1f \n'],mean(temp));
temp=S(HCity);
for ii=1:length(HCity)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii,(StartIndx):EndIndx)));
end

fprintf(' \n \n');

% Aden
temp=Mta(24,(StartIndx):EndIndx);
f=find(temp>0); % Find the weeks where there is at least one attack
fprintf(['Frequency of at least one attack from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Aden: %3.1f \n'],length(temp)./length(f));
fprintf(['Avg. number of attacks from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) 'in Aden: %3.1f \n'],mean(temp));
temp=S(fA==1);
for ii=1:sum(fA)
    fprintf(['Average attack rate per 10,000 from ' datestr(startDateofSim) ' to ' datestr(endDateofSim) ' in ' temp(ii).ADM2_EN ': %5.1f \n'], mean(WIa(ii+3+sum(fS),(StartIndx):EndIndx)));
end