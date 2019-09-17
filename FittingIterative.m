function FittingIterative(XU,tau,AF,CF,RF,dW)
% Runs the fitting for the specified criteria and saves files to folders
% for what is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XU - Specify what you want to use in in the regression model (1 = included, 0 =
%excluded) (1X9)
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates
% tau -Specify a lag of all the factors that are integrated in the model
% (1X6)
    % tau(1) - Past incidence
    % tau(2) - Product of incidence and attacks
    % tau(3) - Product of incidence and conflict
    % tau(4) - Product of incidence and rainfall
    % tau(5) - Perciptiation only
    % tau(6) - Residual incidence
% AF -Specify the attack function to be used
        % AF=0 attack only has effect before;
        %AF=1 Attack has effect only after; 
        %AF=2; Attack has effect before and after
% CF - Specify the conflict function to be used
    % CF=0 linear effect; 
    % CF=1 Hill function with n=1;
    % CF=2; Full hill function; 
% RF- Specify the rainfall function to be used
    % RF=0 Increased incidence for low-rainfall; 
    % RF=1 increased incidence for high rainfall;
    % RF=2 increased incidence for high and low rain fall
% dW- The number of Weeks to be added to the trainign dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Run alorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the data
load('Yemen_Gov_Incidence.mat'); % Incidence data
VS=2.*dW;
RR=length(IData(:,1))-VS;
FWmin=16;
NI=floor((RR-FWmin)./dW)+1; % We add one as the fitting will starts at the minimum FWmin=16. So if there we had to go to index 17 and dW=1, we need two step-sizes since the index will start at 1 for the loop
AIC=zeros(NI,1); % AIC score for the fitting
RPro=zeros(NI,1); % R^2 value for the validation
SSEP=zeros(NI,1); % Sum Square Error for projection
AICp=zeros(NI,1); % AIC for the projection
S = shaperead([ pwd '\ShapeFile\yem_admbnda_adm1_govyem_mola_20181102.shp']); % Shape file for Yemen
load('Conflict_Yemen_Time_Location_Project_Forward.mat'); %Conflict data
WI=IData'; % Transpose the data set such that the number of areas is the row
Ctv=GLevelConflict(ProC,S,length(WI(1,:))); % Send conflict data (time, latitude, longatude) and shapefile of area wanted to catagorize
%% NEED TO ENSURE THAT THE ROWS ARE THE PROPER GOVERNORATES %%
load('Precipitation_Gov_Project_Forward.mat'); % the perceipatiatino data for the differetn areas
Rtv=Rtv(:,1:length(WI(1,:))); % truncate the railfall data to the time period spceified
load('Attack_Yemen_Time_Location_Project_Forward.mat');
tA=GLevelConflict(ProA,S,length(WI(1,:)));% % Send attack data (time, latitude, longatude) and shapefile of area wanted to catagorize
tA(tA>1)=1; % Find the weeks where there are more than one attack that occurs and assume that there is only a single attack

% Load population density
load('PopulationDensity_Yemen.mat'); % loads the population size of the govenerorates (Socotra as the citation did not have their numbers)
P=P';
% Load Health facility density
H=zeros(size(P));

%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
%% Adjust ascpects of functions and data for the fitting

maxtau=4; % The maximum lag allowed for the model
tau(tau>maxtau)=maxtau; % Set to ensure that the lag does not exceed the maximum lag specified
f=find(XU(4:end)==0); % Find the variables not included so the lag can be set to zero
tau(f)=0; % set the lag for components not included to zero

% Specify the lower bounds for the estimated parameters
% beta(1:9)=XU.*[x(1) 10.^(x(2:9))]
% 
% %Attack associated paramters
% DB=10.^x(10);
% DA=10.^x(11);
% % Conflict associated paramters
% K=10.^x(12);
% n=10.^x(13);
% % Rainfall assocaited paramters
% rl=10.^x(14);
% rh=10.^x(15);


lb=[-inf.*ones(1,15)]; % ensuring the lower bound is zero for the last nine paramters and do a log-10 transform to improve searching of paramter space
ub=log10([inf inf inf inf inf inf inf inf inf 1 1 inf inf 10 10]); % specify the upperbound for the parameters 


%% Run the fitting algorithm
NS=1000; % The number of initial starting points used in the fitting process

if isfile([pwd '\Tables\ModelFitIterativeAIC.dat'])
     fmf = fopen([pwd '\Tables\ModelFitIterativeAIC.dat'],'a'); % Writing to the file name
     fmf2 = fopen([pwd '\Tables\ModelFitIterativeAICp.dat'],'a'); % Writing to the file name
     fmf3 = fopen([pwd '\Tables\ModelFitIterativeSSEp.dat'],'a'); % Writing to the file name
     fmf4 = fopen([pwd '\Tables\ModelFitIterativeR2p.dat'],'a'); % Writing to the file name
else
     fmf = fopen([pwd '\Tables\ModelFitIterativeAIC.dat'],'w'); % Writing to the file name
     fmf2 = fopen([pwd '\Tables\ModelFitIterativeAICp.dat'],'w'); % Writing to the file name
     fmf3 = fopen([pwd '\Tables\ModelFitIterativeSSEp.dat'],'w'); % Writing to the file name
     fmf4 = fopen([pwd '\Tables\ModelFitIterativeR2p.dat'],'w'); % Writing to the file name
     fprintf(fmf, 'Model \t Lag_Incidence \t Lag_Attack \t Lag_Conflict \t Lag_RainFall \t Lag_Precipitation \t Lag_External \t AttackFunction \t ConflictFunction \t RainfallFunction '); % The header for the table
     fprintf(fmf2, 'Model \t Lag_Incidence \t Lag_Attack \t Lag_Conflict \t Lag_RainFall \t Lag_Precipitation \t Lag_External \t AttackFunction \t ConflictFunction \t RainfallFunction '); % The header for the table
     fprintf(fmf3, 'Model \t Lag_Incidence \t Lag_Attack \t Lag_Conflict \t Lag_RainFall \t Lag_Precipitation \t Lag_External \t AttackFunction \t ConflictFunction \t RainfallFunction '); % The header for the table
     fprintf(fmf4, 'Model \t Lag_Incidence \t Lag_Attack \t Lag_Conflict \t Lag_RainFall \t Lag_Precipitation \t Lag_External \t AttackFunction \t ConflictFunction \t RainfallFunction '); % The header for the table
     for ii=1:NI
        fprintf(fmf,'\t Iter_%d ',ii);         
        fprintf(fmf2,'\t Iter_%d ',ii);         
        fprintf(fmf3,'\t Iter_%d ',ii);         
        fprintf(fmf4,'\t Iter_%d ',ii);         
     end
     fprintf(fmf,'\n');
     fprintf(fmf2,'\n');
     fprintf(fmf3,'\n');
     fprintf(fmf4,'\n');
end
for ii=1:NI
    FW=FWmin+dW*(ii-1);
    
    opts= optimset('MaxIter',10^6,'MaxFunEvals',10^6,'TolFun',10^(-12),'TolX',10^(-12),'UseParallel',true,'Display','off');  % Specifies the conditions that will be used in the fitting processs
    problem = createOptimProblem('lsqnonlin','objective',@(x)OFunc(x,WI(GNZI,1:FW),tA(GNZI,1:FW),Ctv(GNZI,1:FW),Rtv(GNZI,1:FW),XU,tau,maxtau,AF,CF,RF,P(GNZI),H(GNZI)),'x0',log10(rand(15,1)),'lb',lb,'ub',ub,'Aineq',[],'bineq',[],'Aeq',[],'beq',[],'options',opts);
    ms = MultiStart('UseParallel','always'); % specifies that we run the algorithm in parallel
    [par,fval] = run(ms,problem,NS); %starts at NS random initial points to thoroghly search the paramter space
    % Evaluate the number of paramters that are being used in the estimation 
    [k,beta,DB,DA,K,n,rl,rh]=RetParameter(par,XU,AF,CF,RF);

    % Calculate the AIC score for the model fit
    nData=FW-maxtau; % The number of data points available for the model fitting
    AIC(ii)= AICScore(k,nData,fval); % Calculate the AIC score
        


    %% Set up information for the projection

    %% Run the projection
    [Yt,Pt]= ModelProjection(beta,WI(GNZI,1:(FW)),tA(GNZI,1:(FW+VS)),DB,DA,Ctv(GNZI,1:(FW+VS)),K,n,Rtv(GNZI,1:(FW+VS)),RF,rl,rh,tau,maxtau,CF,P(GNZI),H(GNZI),VS); % Run the projection
    %% Compute the R^2 value of the model projection
    PWI=WI(GNZI,(FW+1):(FW+VS));
    RPro(ii)=corr(Pt(:),PWI(:)).^2; % The R^2 value for the projection component of the model fit
    SSEP(ii)=sum((Pt(:)-PWI(:)).^2);
    nData=length(Pt(:)); % The number of data points available for the model fitting
    AICp(ii)= AICScore(k,nData,SSEP(ii)); % Calculate the AIC score  

end
RFF=struct('N',{'None','Low','High','LowHigh'});
CFF=struct('N',{'None','Linear','Hill(n=1)','Hill'});
AFF=struct('N',{'None','Before','After','BeforeAfter'});

fprintf(fmf,['%d  \t %d  \t %d  \t %d  \t %d  \t %d  \t %d  \t ' AFF((AF+1).*XU(5)+1).N '\t ' CFF((CF+1).*XU(6)+1).N ' \t ' RFF((RF+1).*XU(7)+1).N ],[XU*((2.^[0:(length(XU)-1)])') tau.*XU(4:end)] );
for ii=1:NI   
    fprintf(fmf,' \t %6.2f',AIC(ii)); 
end
fprintf(fmf,'\n');


fprintf(fmf2,['%d  \t %d  \t %d  \t %d  \t %d  \t %d  \t %d  \t ' AFF((AF+1).*XU(5)+1).N '\t ' CFF((CF+1).*XU(6)+1).N ' \t ' RFF((RF+1).*XU(7)+1).N ],[XU*((2.^[0:(length(XU)-1)])') tau.*XU(4:end)] );
for ii=1:NI   
    fprintf(fmf2,' \t %6.2f',AICp(ii)); 
end
fprintf(fmf2,'\n');


fprintf(fmf3,['%d  \t %d  \t %d  \t %d  \t %d  \t %d  \t %d  \t ' AFF((AF+1).*XU(5)+1).N '\t ' CFF((CF+1).*XU(6)+1).N ' \t ' RFF((RF+1).*XU(7)+1).N ],[XU*((2.^[0:(length(XU)-1)])') tau.*XU(4:end)] );
for ii=1:NI   
    fprintf(fmf3,' \t %6.2f',SSEP(ii)); 
end
fprintf(fmf3,'\n');


fprintf(fmf4,['%d  \t %d  \t %d  \t %d  \t %d  \t %d  \t %d  \t ' AFF((AF+1).*XU(5)+1).N '\t ' CFF((CF+1).*XU(6)+1).N ' \t ' RFF((RF+1).*XU(7)+1).N ],[XU*((2.^[0:(length(XU)-1)])') tau.*XU(4:end)] );
for ii=1:NI   
    fprintf(fmf4,' \t %4.3f',RPro(ii)); 
end
fprintf(fmf4,'\n');

fclose('all');
end