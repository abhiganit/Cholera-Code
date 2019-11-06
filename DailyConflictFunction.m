function F= DailyConflictFunction(DateV,C,Func,X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DateV - Vecotr of the date interested in the output
% C - 1 x N vector containing index of day event occured in the area of
% interest starting Sept. 30 2016 (i.e. index 1 is 
% Func - The function to evaluate the conflict
% X the threshold that will be used in the calculation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F - the value of the function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comput function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startDateofData = datenum('09-30-2016');% Start date of the daily conflict data
IndX=DateV-startDateofData; % Subtract startDateofData to obtain the index and subtracty one. Since if we wanted Oct 1, 2016 it would be any with index 1
if(Func>0)
    Ct=C-X; % Adjust the treshhold
else
 Ct=C;
end
Y=Ct(:,IndX);
Y(Y<0)=0;
if(Func==2)
    Y(Y>1)=1;
end

F=mean(Y,2);

end

