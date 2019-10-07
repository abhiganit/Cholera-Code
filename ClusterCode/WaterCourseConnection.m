function Mt = WaterCourseConnection(W,S)

Temp=zeros(length(W),length(S));

for ii=1:length(S)
    parfor jj=1:length(W)
      in=inpolygon(W(jj).X,W(jj).Y,S(ii).X,S(ii).Y);
      f=find(in>0);
      Temp(jj,ii)=min([length(f) 1]);
    end
end

Mt=zeros(length(S),length(S));
for ii=1:length(S)
   f=find(Temp(:,ii)==1) ;
   for jj=1:length(f)
      if(sum(Temp(f(jj),:))>1)
         g=find(Temp(f(jj),:)==1);
         h=find(g~=ii);
         h=g(h);
         for zz=1:length(h)
            Mt(ii,h(zz))=1+Mt(ii,h(zz));
            Mt(h(zz),ii)=1+ Mt(h(zz),ii);
         end
      end
   end
end
end