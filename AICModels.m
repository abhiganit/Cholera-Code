load('Yemen_Gov_Incidence.mat'); % Incidence data
PDS=0.8;
WI=IData'; % Transpose the data set such that the number of areas is the row
maxtau=4;
%Find areas where we have non-zero incidence over course of epidemic
GNZI=find(sum(WI,2)~=0); % Critical if we are estimating beta_0 otherwise does not make a difference
NWF=153;
NGS=floor(length(GNZI)*PDS);
Itemp=sum(WI(GNZI,:),2);
GTF=zeros(length(NGS),1); % We use the top and bottom gov wrt incidence in the fitting of the model and cross validate to the remaining ones in the middle
% Find the top max
for ii=1:ceil(NGS/2)
f=find(Itemp==max(Itemp)); % Find the maximum
Itemp(f)=0; % set maximum to zero for it is no longer selected
GTF(ii)=f; % Record index
end
Itemp=sum(WI(GNZI,:),2); % Recalc number of cases
% Find the minimum contributors
for ii=(ceil(NGS/2)+1):NGS
f=find(Itemp==min(Itemp)); % Select minimum
Itemp(f)=max(Itemp); % Set to maximum for not selected again
GTF(ii)=f; % Record index
end
GTF=sort(GTF)';
WI=WI(GNZI(GTF),(1+maxtau):NWF);
N=length(WI(:));

%% Compute AIC for the models starting from the Incidence model

AIC=zeros(4,1);

load('ForwardSelectionNoConflictNoRain-PercentDataSet=80-alpha=0.mat')
AIC(1)=AICScore(kv(end),N,RSSv(end));


load('ForwardSelectionNoRain-PercentDataSet=80-alpha=0.mat')
AIC(2)=AICScore(kv(end),N,RSSv(end));


load('ForwardSelectionNoConflict-PercentDataSet=80-alpha=0.mat')
AIC(3)=AICScore(kv(end),N,RSSv(end));

load('ForwardSelection-PercentDataSet=80-alpha=0.mat')
AIC(4)=AICScore(kv(end),N,RSSv(end));

dAIC=AIC-min(AIC)