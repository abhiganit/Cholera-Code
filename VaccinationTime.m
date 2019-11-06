function [V1,V2]= VaccinationTime(Scale,NW)
% Returns matrix of the number of doses given over the course of the
% outbreak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale - Scale=1 indicates the governorate level and scale=2 is the
% district level
% NW - the number of weeks to be included in the matrix relative to Oct 3,
% 2016 (i.e. NW =1 then only week of Oct 3 2016 is included)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1 - Matrix of people to receive at least one dose
% V2 - Matrix of people to receive at least two doses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startDateofSim = datenum('10-03-2016');% Start date
if(Scale==1)
    S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
    S={S.ADM1_EN};
    V1=zeros(22,NW);   
    V2=zeros(22,NW);   
    
    % Vaccination in Aden
    VaxDate=datenum('05-06-2018'); % May 6 2018
    indx=ceil((1+VaxDate-startDateofSim)./7);
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Aden'},S(ii)); 
    end
    V1(ff==1,indx)=306462;    
    V2(ff==1,indx)=235833;
    
    
    % Vaccination in Al Hodeidah and Ibb
    
    VaxDate=datenum('08-04-2018'); % August 4 2018 
    indx=ceil((1+VaxDate-startDateofSim)./7);
    % Al Hodeidah
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Hudaydah'},S(ii)); 
    end
    V1(ff==1,indx)=312652;
    V2(ff==1,indx)=185026;
    
    % Ibb
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Ibb'},S(ii)); 
    end
    V1(ff==1,indx)=105886;
    V2(ff==1,indx)=93429;
    
    % Vaccination Taizz, Al Dhale and Aden
    
    VaxDate=datenum('02-24-2019'); % Feb 24, 2019
    indx=ceil((1+VaxDate-startDateofSim)./7);
    
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Taizz'},S(ii)); 
    end
    V1(ff==1,indx)=137803;
    
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Dhale''e'},S(ii)); 
    end
    V1(ff==1,indx)=254462;
    
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Aden'},S(ii)); 
    end
    V1(ff==1,indx)=86270;
    
    % Vaccination Amana Al Asmiah
    
    VaxDate=datenum('04-24-2019'); % Feb 24, 2019
    indx=ceil((1+VaxDate-startDateofSim)./7);
    
    ff=zeros(length(S),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Amanat Al Asimah'},S(ii)); 
    end
    V1(ff==1,indx)=1088018;
elseif(Scale==2)    
    
    
    V1=zeros(24,NW);   
    V2=zeros(24,NW);   
    SD = shaperead([ pwd '\ShapeFile\yem_admbnda_adm2_govyem_mola_20181102.shp']); % Shape file for Yemen
    
    
    fS=zeros(length(SD),1);
    for ii=1:length(fS)
      fS(ii)=strcmp({'Amanat Al Asimah'},SD(ii).ADM1_EN); 
    end
    fS=find(fS==1);

    fA=zeros(length(SD),1);
    for ii=1:length(fA)
      fA(ii)=strcmp({'Aden'},SD(ii).ADM1_EN); 
    end

    fA=find(fA==1);


    SD={SD([29 31 71 fS' fA']).ADM2_EN};
    
    % Vaccination in Aden
    VaxDate=datenum('05-06-2018'); % May 6 2018
    indx=ceil((1+VaxDate-startDateofSim)./7);
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Buraiqeh'},SD(ii)); 
    end
    V1(ff==1,indx)=79447;    
    V2(ff==1,indx)=49051;
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Attawahi'},SD(ii)); 
    end
    V1(ff==1,indx)=59589;    
    V2(ff==1,indx)=43717;
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Mualla'},SD(ii)); 
    end
    V1(ff==1,indx)=50534;    
    V2(ff==1,indx)=44556;
    
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Khur Maksar'},SD(ii)); 
    end
    V1(ff==1,indx)=52960;    
    V2(ff==1,indx)=38316;
    
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Craiter'},SD(ii)); 
    end
    V1(ff==1,indx)=63932;    
    V2(ff==1,indx)=60193;
    
    % Vaccination in Al Hodeidah
    
    VaxDate=datenum('08-04-2018'); % August 4 2018 
    indx=ceil((1+VaxDate-startDateofSim)./7);
    % Al Hodeidah
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Hali'},SD(ii)); 
    end
    V1(ff==1,indx)=144991;
    V2(ff==1,indx)=49410;
    
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Marawi''ah'},SD(ii)); 
    end
    V1(ff==1,indx)=167661;
    V2(ff==1,indx)=135616;
    
    %  Aden
    
    VaxDate=datenum('02-24-2019'); % Feb 24, 2019
    indx=ceil((1+VaxDate-startDateofSim)./7);
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Dar Sad'},SD(ii)); 
    end
    V1(ff==1,indx)=86270;
    
    % Vaccination Amana Al Asmiah
    
    VaxDate=datenum('04-24-2019'); % Feb 24, 2019
    indx=ceil((1+VaxDate-startDateofSim)./7);
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Shu''aub'},SD(ii)); 
    end
    V1(ff==1,indx)=293758;
    
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'As Sabain'},SD(ii)); 
    end
    V1(ff==1,indx)=660652;
    
    
    ff=zeros(length(SD),1);
    for ii=1:length(ff)
        ff(ii)=strcmp({'Al Wahdah'},SD(ii)); 
    end
    V1(ff==1,indx)=133608;
    
    % Hodeidah City
    V1(22,:)=sum(V1(1:3,:),1);
    V2(22,:)=sum(V2(1:3,:),1);
    
    
    % Sana'a City
    V1(23,:)=sum(V1(3+[1:length(fS)],:),1);
    V2(23,:)=sum(V2(3+[1:length(fS)],:),1);
    
    
    % Aden
    V1(24,:)=sum(V1(3+length(fS)+[1:length(fA)],:),1);
    V2(24,:)=sum(V2(3+length(fS)+[1:length(fA)],:),1);
end

end

