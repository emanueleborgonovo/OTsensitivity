function [Popt,fopt]=transsimp2(a,b,C)
% Transportation simplex as of Luenberger Yu, Jarre Stoer
% [P,COST]=TRANSSIMP2(A,B,C) computes the optimal transport for
% the cost matric C with marginals A and B
%
% next try with tree structure
[n,m]=size(C); % cmax=max(C(:));
verbose= false; % true; % 
initfeas='NW2';  % 'Vogel'; % 'NW'; % 
assert(size(a,1)==n && size(b,1)==m,'mass and cost matrix dimension mismatch')
%assert(abs(sum(a)-sum(b))<max(n,m)*eps(cmax), 'No feasible solution for the transport problem.');

% preallocate working arrays
nm1=n+m-1;
iB=zeros(nm1,2); % indices of basis
sP=zeros(nm1,1);   % sparse solution
sC=zeros(nm1,1);   % sparse costs
u=nan(n,1);v=nan(m,1); % simplex multipliers
tree=zeros(nm1,1); 
wayBack=zeros(nm1,1); 
orphans=false(nm1,1);
visited=false(nm1,1);

%% Step 1. Compute an initial basic feasible solution using the Northwest Corner
% Rule or some other method.
switch initfeas
case 'NW'
%
 i=1;j=1;
 supply=a(i);demand=b(j);
 tree(1)=0; % root
 for cnt=1:(n+m-1)
  iB(cnt,:)=[i,j];
  sC(cnt)=C(i,j);
  if(supply>demand) 
   sP(cnt)=demand; 
   supply=supply-demand;
   if(j<m) 
   % second coordinate changes, first one stays the same
    j=j+1;  
    tree(cnt+1)=cnt; 
   end
   demand=b(j); 
  else 
   sP(cnt)=supply; 
   demand=demand-supply;
   if(i<n)
    i=i+1; 
    tree(cnt+1)=-cnt;
   end
   supply=a(i); 
  end
 end
 root=1; 
case 'NW2'
 % Northwestern Rule, alternative
 i=1;j=1;
 supply=b(j);demand=a(i);
 tree(1)=0; % root
 for cnt=1:(n+m-1)
  iB(cnt,:)=[i,j];
  sC(cnt)=C(i,j);
  if(supply>demand) 
   sP(cnt)=demand; 
   supply=supply-demand;
   if(i<n)
    % first coordinate changes
    i=i+1; 
    % second coordinate stays the same
    tree(cnt+1)=-cnt; 
   end
   demand=a(i); 
  else 
   sP(cnt)=supply; 
   demand=demand-supply;
   if(j<m)
    j=j+1; 
    tree(cnt+1)=cnt;
   end
   supply=b(j); 
  end
 end
 root=1;
 case 'Vogel'
%% Vogel Approximation Method, second try
 as=a;bs=b;%c=C;
 bcnt=0; % # base entries
 rowsUsed=0;
 colsUsed=0;
 n0=(1:n)';
 m0=(1:m)';

 while rowsUsed<n-1 && colsUsed<m-1 % (any(as))
     rowsUsed+colsUsed
  c=C(n0,m0);
  [rmax,l]=max(diff(mink(c,2,1)));   % minimal and minimal-but-one
  [cmax,k]=max(diff(mink(c,2,2)')); 
  if(rmax>=cmax)
   [~,k]=min(c(:,l));
  else
   [~,l]=min(c(k,:));
  end
  i=n0(k);j=m0(l);
  
  aa=as(i);bb=bs(j);
  d=aa-bb;
  if d>0
   as(i) = d;
   bs(j) = 0;
 
    bcnt=bcnt+1;
    sP(bcnt)=bb; 
    iB(bcnt,:)=[i,j];    

    m0(m0==j)= [];
    colsUsed=colsUsed+1;
  else %if d<0
    as(i) = 0;
    bs(j) =-d;

    bcnt=bcnt+1;
    sP(bcnt)=aa;
    iB(bcnt,:)=[i,j];    

    n0(n0==i)= [];
    rowsUsed=rowsUsed+1;
  end
 end
% Dempe/Schreier: Existiert nur noch eine ungestrichene Zeile oder eine ungestrichene
% Spalte, dann gehe zur Nordwesteckenregel über.
 % either n0 or m0 will be a scalar
 iB(bcnt+1:end,:)=[n0*ones(size(m0)),m0*ones(size(n0))];
 error('no tree implemented.')
otherwise
 error('Unknown initialization method');   
end % select method

while (true)
%% Step 2. Compute the simplex multipliers and the relative cost coefficients. 
  % Simplex multipliers
  % ********************************************* 
  % Step 1. Assign an arbitrary value tor any one of the multipliers.

% Jarre / Stoer 
% Wir können ohne Einschränkung u1 = 0 wählen.
% Für die Nachbarn Dl von S1, die mit S1 durch eine Kante in G(J ) verbunden sind, folgen
% aus u1 = 0 und u1 +v ell = c1ell die Werte vell . Für deren Nachbarn Sk (in G(J )) sind dann
% wiederum die Werte uk durch uk +vell = ckell eindeutig gegeben. So lassen sich sukzessive
% alle ui und v j bestimmen.

   % lastleaf and tree are initialized
  visited(:)=false; %(nm1,1);
   
   % built a tree traverser (depth first) so that 
   % it is assured that the parent nodes are already visited
   % the traversal is implemented by an indirect index lookup 
   % for the for loop, which makes the endfor spending 5% of total time
  visited(root)=true;
   % find the leaves
  leavs=0:nm1;
  leavs(abs(tree)+1)=[];
   % traverse back until a visited element is hit
  wayBack(1)=root;
  ii=1:nm1; % index reverser
  lower=2;upper=2;
  for l=leavs
   k=l;
   while ~visited(k)
    visited(k)=true;
    wayBack(upper)=k;
    upper=upper+1;k=abs(tree(k));
   end
   % each segment, starting from a leaf and ending at an element with 
   % parent node root or previously visited, is traversed in reverse order
   if(lower+1<upper)
    ii(lower:upper-1)=upper-1:-1:lower; % revert
   end
   lower=upper;
  end
  
  for k=wayBack(ii)' % compute the multiplieres
   q=tree(k);
   % horizontal or vertical link?
   if(q>0)
    v(iB(k,2))=sC(k)-u(iB(q,1));
   elseif(q<0)
    u(iB(k,1))=sC(k)-v(iB(-q,2));
   else%(q==0)
    % root node, first entry
    u(iB(k,1))=0; % we might reuse the old value here and update only the changed part of the tree
    v(iB(k,2))=sC(k);
   end
  end
  if(verbose)
   disp('Simplex multipliers');
   u',v'
  end

  % **********************************************
% Profiler: 8% spent here
  R=C-u-v'; % relative costs
% Profiler: 4% spent here -- remove?
  % enforce exact zeros for all indices in iB
  for i=1:nm1
   R(iB(i,1),iB(i,2))=0; 
  end
%% Step 3. Select a nonbasic variable corresponding to a negative cost coefficient to
% enter the basis (usually the one corresponding to the most negative cost coefficient).
  [cr,ij]=min(R(:)); % flatten for speed
% If all relative cost coefficients are nonnegative, stop; the solution is optimal. 
  if(cr>-10*eps), break; end
  % selected variable index (don't use i,j as iterators, though)
  i=1+mod(ij-1,n);j=1+(ij-i)/n; 

 % Compute the cycle of change
  bigcycle=[];
  cycle=[];
% first match suffices, as cycle is cleaned later 
  matchI=find(iB(:,1)==i,1); 
  matchJ=find(iB(:,2)==j,1);

  % reroot at horizontal match
  k0=matchI;
  q=tree(k0);
  if(q~=0) % something to do
   tree(k0)=0;
   while(q~=0)
     k=abs(q);s=sign(q);
     q=tree(k); tree(k)=s*k0;k0=k;
   end 
  end
  % traverse back from vertical match until we hit the (new) root
  k=matchJ;
  % also record the inverse tree
  eert=tree;
  bigcycle(1)=k;
  q=eert(k); eert(k)=0; k0=k;
  while q~=0
   k=abs(q);s=sign(q);
   bigcycle(end+1)=k;
   q=tree(k); eert(k)=s*k0;k0=k;
  end

  % skip over dead-end leaf for small cycle
  % ...
  B=[iB(bigcycle,:);i,j];
  horv=2; % horizontal or vertical, start with vertical
  uv=B(1,horv);
  cycleindx=[];
%  thtanew=inf;
  for l=1:length(bigcycle)

    if B(l+1,horv)~=uv
	 cycle(end+1)=bigcycle(l);
	 cycleindx(end+1)=l;
	 horv=3-horv; % 1 or 2
	 uv=B(l,horv);
	 % perform search for thta already here?
%     if(horv==1)
%         if(uv<thtanew), thtanew=uv; lnew=length(cycleindx);end
%     end
    end
  end   
  if(verbose)
   disp('Cycle index (new coord last)');
   disp(B') % the cycle index
  end

  % cycle and tree wind differently 
 assert(iB(cycle(end),1)==i && iB(cycle(1),2)==j,'Pivot cycle winds wrong way.')  
%% Set theta equal to the smallest basic variable
% with a minus assigned to it. 
  [thta,l]=min(sP(cycle(1:2:end)));
  l=1+2*(l-1); % back to original coordinates
  k=cycle(l);  
  
  % update tree structure
  uncycle=true(nm1,1);
  uncycle(bigcycle)=false;
  
  orphans(:)=false;
  if(l==1) % remove first element in the cycle (matchJ, if not cleaned)
    % pointing to k with first coordinate?
   orphans(uncycle)=iB(uncycle,1)==iB(k,1);
  
   if(any(orphans)), tree(orphans)=bigcycle(1+cycleindx(l)); end % or bigcycle
   kk=k;
   while kk~=0 % matchJ
      q=eert(kk);tree(kk)=q;kk=abs(q);
   end
   
   if(k~=matchJ) % chopped off from bigcycle
    tree(matchJ)=-k;
   end
   tree(matchI)=k; %matchI was current root
   root=k; %matchJ;
   tree(root)=0;

  elseif(l==length(cycle))  % remove last element in the cycle (matchI)
  % pointing to k with first coordinate?
   orphans(uncycle)=iB(uncycle,2)==iB(k,2);
   if(any(orphans)), tree(orphans)=-bigcycle(cycleindx(l)-1); end
   tree(k)=-matchJ;
   if(k~=matchI) % chopped off from bigcycle
    tree(matchI)=cycle(end);
   end
   root=bigcycle(cycleindx(l)-1); % eert(k); % cycle(l-1); 
   tree(root)=0;
  else % cycle(l-1) and cycle(l+1) exist
   orphans(uncycle)=iB(uncycle,1)==iB(k,1); % i;
   if(any(orphans))
   if(verbose), fprintf('horizontal orphans: %d\n',find(orphans)); end
   tree(orphans)=cycle(l+1); end % bigcycle(1+cycleindx(l)); 
   orphans(:)=false;
   orphans(uncycle)=iB(uncycle,2)==iB(k,2); % j;
   if(any(orphans))
      if(verbose), fprintf('vertical orphans: %d\n',find(orphans)); end
   tree(orphans)=-cycle(l-1);  end %-bigcycle(cycleindx(l)-1); 
   kk=k;
   while kk~=0
      q=eert(kk);tree(kk)=q;kk=abs(q);
   %if(kk==0), break; end % matchJ
   end
   % chopped off from bigcycle
   if(cycle(1)~=matchJ)
    tree(cycle(1))=-matchJ;
   end
   if(cycle(end)~=matchI) 
    tree(matchI)=cycle(end);
   end
   tree(k)=-matchJ;
   tree(matchI)=k;
   root=matchJ; 
  end

if verbose
  disp('Updated: pivot,val');
k
thta
end
%Update the solution.
if(thta~=0.0), sP(cycle)=sP(cycle)+thta*(-1).^(1:length(cycle))'; end
sP(k)=thta;
sC(k)=C(i,j);
iB(k,:)=[i,j]; % update basis

if(verbose)
 disp('cost, potential');
 disp([sC,sP]');
end
end % while
% construct full matrix
Popt=zeros(n,m);
for i=1:(n+m-1)
Popt(iB(i,1),iB(i,2))=sP(i);
end
fopt=u'*a+v'*b;
end % function

function xxtest
 b=[1,5,2,8,2]';
 a=[3,8,1,6]';
 C=[3 4 6 8 9; 2 2 4 5 5; 2 2 2 3 2; 3 3 2 4 2];
[p,q]=transsimp2(a,b,C)
 transsimp2(b,a,C')
 %% always quadratic
 intinf=9; % integer infinity (unattractive cost)
 D=[intinf*ones(length(b)),C';C,intinf*ones(length(a))];
 d=[b;a];
 [q,p]=transsimp2(d,d,D);q,p/2
 %%
end