function RC = RebelControl(R,S)
% Determines whether any part of the polygon in S is within the boundary of
% Rebel control 

RC=zeros(length(S),1);

for rr=1:length(RC)
   in=inpolygon(S(rr).X,S(rr).Y,[R.X],[R.Y]);    
   in2=inpolygon([R.X],[R.Y],S(rr).X,S(rr).Y);  
   inc=min(max(in)+max(in2),1);
   RC(rr)=max(inc);
end
end

