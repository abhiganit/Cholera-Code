function XS = RetStartpoint(XU,tau,AF,CF,RIF,RF,lb,ub)
% RetStartpoint(XU,tau,AF,CF,RIF,RF) returns a suitable starting point for
% the algorithm based on some of the past model fits
%   Detailed explanation goes here
% XU - Specify what you want to use in in the regression model (1 = included, 0 =
%excluded) (1X11)
        % XU(1)- beta_0
        % XU(2) - population density
        % XU(3) - number of health facilities 
        % XU(4) - Past incidence
        % XU(5) - Product of incidence and attacks
        % XU(6) - Product of incidence and conflict
        % XU(7) - Product of incidence and rainfall
        % XU(8) - Rainfall only        
        % XU(9) - Incidence in other govnorates
        % XU(10)- Attacks only
% tau -Specify a lag of all the factors that are integrated in the model
% (1X10)
    % tau(1) - population density incidence
    % tau(2) - health zone incidence
    % tau(3) - Past incidence
    % tau(4) - Product of incidence and attacks
    % tau(5) - Product of incidence and conflict
    % tau(6) - Product of incidence and rainfall
    % tau(7) - Perciptiation only
    % tau(8) - Incidence in other govneroates
    % tau(9)- Attack only
    % tau(10)- Rebel control
% AF -Specify the attack function to be used
        % AF=0 attack only has effect before;
        %AF=1 Attack has effect only after; 
        %AF=2; Attack has effect before and after
% CF - Specify the conflict function to be used
    % CF=0 linear effect; 
    % CF=1 Hill function with n=1;
    % CF=2; Full hill function; 
% RIF- Specify the rainfall function to be used for rainfall*incidence
% covariate
    % RIF=0 Increased incidence for low-rainfall; 
    % RIF=1 increased incidence for high rainfall;
    % RIF=2 increased incidence for high and low rain fall
% RF- Specify the rainfall function to be used for rainfall covariate
    % RF=0 Increased incidence for low-rainfall; 
    % RF=1 increased incidence for high rainfall;
    % RF=2 increased incidence for high and low rain fall
% lb- the log_10 lower bounds
% ub- the log_10 upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% XS- The starting points 
% Specify the lower bounds for the estimated parameters
% beta=[10.^x(1:11)'].*XU;
% 
% %Attack associated paramters
% DB=10.^x(11);
% DA=10.^x(12);
% % Conflict associated paramters
% K=10.^x(13);
% n=10.^x(14);
% % Rainfall assocaited paramters
% rl=10.^x(15);
% rh=10.^x(16);
% DAE=10.^x(17);

% Read table of past fits
F=struct('N',{'None','Low','High','LowHigh'});
CFF=struct('N',{'None','Linear','Linear_Threshold','Poly_Threshold'});
AFF=struct('N',{'None','Before','After','BeforeAfter'});

T=readtable([pwd '\Tables\ProjectionModelFitSummary.dat']);
if(~isempty(T))
    XUt=repmat([XU(1) XU(2:end).*tau],length(T.Intercept),1);
    XT=[T.Intercept    T.Population    T.HealthFacilities    T.Incidence    T.Attack_Incidence    T.Conflcit_Incidence    T.Rainfall_Incidence    T.Precipitation    T.External_Incidence    T.Attack T.RebelControl];
    Ob=T.ObjectiveValue;
    dX=XUt-XT;
    mindX=min(dX')';
    f=find(mindX>=0);
    ftemp=find(contains(T.Attack_Function(f),AFF(AF*XU(5)+1).N));
    f=ftemp;
    ftemp=find(contains(T.Conflict_Function(f),CFF(CF*XU(6)+1).N));
    f=ftemp;
    ftemp=find(contains(T.Rainfall_Incidence_Function(f),F(RIF*XU(7)+1).N));
    f=ftemp;
    ftemp=find(contains(T.Precipitation_Function(f),F(RF*XU(8)+1).N));
    f=ftemp;
    if(isempty(f))
        XS=lb+(ub-lb).*rand(size(lb));
    else
        g=find(Ob(f)==min(Ob(f)),1);
        indx=f(g);
        XUs=XT(indx,:);
        taus=XUs(2:end);
        XUs(XUs>1)=1;
        load([pwd '\Tables\ProjectionModel-XU=' num2str(XUs*((2.^[0:(length(XUs)-1)])')) '-CF=' num2str(CF*XUs(6)) '-AF=' num2str(AF*XUs(5)) '-RIF=' num2str(RIF*XUs(7)) '-PF=' num2str(RF*XUs(8)) '-tau=' num2str(taus(XUs(2:end)>0)) '.mat'],'par');
        f=find(XUs==0);
        par(f)=lb(f)+1;
        XS=par';
    end
else    
    XS=lb+(ub-lb).*rand(size(lb));
end
XS=XS';
end

