function [c,bi]=lapmack(C)
% LAPMACK Linear Assignment Problem using Mack's approach (as discussed by
% Jonker Volgenant)

%% DEFUNCT!

n=size(C,1);
dsum=0;
% step 1: Determine the basis in the cost matrix
[ci,bi]=min(C,[],2); % bi: base_i

cnt=accumarray(bi,1,[n,1]);
iter=0;
ba=zeros(n,1); % alternative basis, unpopulated
%forl=1;
while(true)
 S=ones(n,n);for i=[1:n;bi']; S(i(1),i(2))=-1; end
 C.*S
 pause
% step 3: select a column which contains more than one base
 lcol=false(n,1); % logical vector
 %if(forl==1)
 %    firstorlast='first';
 %    forl=0;
 %else
 %    firstorlast='last';
 %    forl=1;
 %end
 col=find(cnt>1,1); %,firstorlast); 
% step 2: Termination if every column contains one base
 if(isempty(col))
   c=-dsum; % cost
   for i=1:n; c=c+C(i,bi(i)); end    
   return 
 end
 lcol(col)=true;
 row=find(bi==col);

 while true

  iter=iter+1; % max. n steps
%  if iter>n, warning('oops'); return; end 
%% step 4
 %  for i=row'
 % set m = min { C_i,k-C_i,base_i | k in col}
 % determine delta=min (mi| i in ROW}
 % Let rr be a row and kk a column for which delta is assumed
  ncol=find(~lcol);
  [m,j]=min(C(row,~lcol)-ci(row),[],2);
  [d,r]=min(m);
  rr=row(r);
  kk=ncol(j(r));
  
 % step 5: Audjust the dual solution
  C(:,lcol)=C(:,lcol)+d;dsum=dsum+d*sum(lcol);
  if cnt(kk)==0 % kk contains no base
   break; % out of step 4 loop
  else
 % step 7: column kk is base for some row
  % mark column kk as alternative base for row rr
   ba(rr)=kk;
  % col= col u {kk} row=row u { i| base_i=kk}
   %col(end+1)=kk; 
   lcol(kk)=true;
   bb=find(bi==kk);
   row(end+(1:length(bb)))=bb;
%   row(randperm(length(row)))=row; % perturb to avoid loop formation
  % goto step 4
  end
 end
  % Step 6: Augmentation
 
 % alter the current set of bases along the alternating path, starting in
 % column kk
 bi(rr)=kk; % enter directly
 bi(ba>0)=ba(ba>0);
 ba(:)=0;
 for i=1:n, ci(i)=C(i,bi(i)); end
 cnt=accumarray(bi,1,[n,1]);
 col=find(cnt>1,1); % prepare step 3
 % goto step 2
end
%
end

function test
%%
C=[3 7 6 6 ; 1 6 8 8 ; 3 0 8 1; 0 7 9 9]
lapmack(C)
lapmack([C,C;C,C]) % and bugged
%%
end