function [assignment,cost]=revsimplex(a,b,costm)
% REVSIMPLEX Revised simplex of Luenberger Ye
% port from the transport.R package

 THRESHOLD=-1e-6;
 [m,n]=size(costm);

 assert(length(a)==m && length(b)==n,'Size mismatch.');

 assignment=zeros(m,n);
 basis=false(m,n);

 maxdim=max(m,n);
 iter=0;
 next_row=1; % matlab has base 1

% stretch space
 basis_byrow=zeros(m,n);
 basis_byrow_top=zeros(m,1);

 basis_bycol=zeros(n,m);
 basis_bycol_top=zeros(n,1);

 u=zeros(m,1);
 v=zeros(n,1);

 is_computed_u=false(m,1);
 is_computed_v=false(n,1);

 maxcycle=2*maxdim; % was m+n
 list=zeros( maxcycle,1);
 is_row=zeros( maxcycle,1);
 circlea=zeros( maxcycle,1);
 circleb=zeros( maxcycle,1);

 candlist=zeros(maxdim,1);
 rem_curr=ones(maxdim,1); % indices
 rem_next_branch=ones(maxdim,1); %indices
 rem_do_rowscan=false(maxdim,1);

% initial feasible solution, northwestern rule 
 % disp('In northwestern') %% DBG
 aleft = a;
 bleft = b;
 aisleft = true( m, 1);
 bisleft = true( n, 1);
 degenisj = false( n, 1);
 degennum = 0;
 numleft = m;
 while numleft > 0 
    for i = 1:m     
      if aisleft(i) 
	     mini = inf;
       for j = 1:n 
        if bisleft(j) 
         if costm(i,j) < mini 
          mini = costm(i,j);
          jstar = j;
         end
        end
       end
       mass = min(aleft(i),bleft(jstar));
       aleft(i) = aleft(i) - mass;
       bleft(jstar) = bleft(jstar) - mass;
       assignment(i,jstar) = mass;
       basis(i,jstar) = true;
       if aleft(i) <eps  % == 0
          aisleft(i) = false;
          numleft = numleft - 1;
       end
       if bleft(jstar) < eps % == 0
          bisleft(jstar) = false;
       end
       if aleft(i) <eps && bleft(jstar) <eps % degenerate case
        if degennum > 0
         for j = 1:n
          if degenisj(j)  && ~bisleft(j)
           degenj = j;
           basis(degeni, degenj) = true;
           break;
          end
         end
        end
        degeni = i;
        degenisj = bisleft;
        degennum=degennum+1;
       end
      end
    end
 end
 clear aleft bleft aisleft bisleft degenisj

% init helpers
% disp('In inithelpers') % DBG
 for i = 1:m
   for j = 1:n
    if basis(i,j)
        % stack lists
      basis_byrow_top(i) = basis_byrow_top(i)+1;
      basis_byrow(i,basis_byrow_top(i)) = j;
      basis_bycol_top(j) = basis_bycol_top(j)+1;      
      basis_bycol(j,basis_bycol_top(j)) = i;
    end
   end
 end
 found=true;
 while found
    iter=iter+1;
    
    % update_transport_rowmostneg
    % found = new_basic_variable_rowmostneg(state);
    % disp('in updatetransport'); % DBG
    top = 1; % actually, list seems to be a stack       
    curr = 1;
    is_computed_u(:) = false;  
    is_computed_v(:) = false;
    u(1) = 0.0; is_computed_u(1) = true;
    list(top) = 1;
    is_row(top) = true;
    %top=top+1;
    
    while(curr <= top)
     pos = list(curr);
     if is_row(curr)
      for j = 1:basis_byrow_top(pos)
       newind = basis_byrow(pos, j);
       if ~is_computed_v(newind) 
          v(newind) = costm( pos, newind) - u(pos);
          is_computed_v(newind) = true;
          top = top + 1;
          list(top) = newind;
          is_row(top) = false;
       end
      end
     else 
      for i = 1:basis_bycol_top(pos)
        newind = basis_bycol( pos, i);
        if ~is_computed_u(newind)
          u(newind) = costm( newind, pos) - v(pos);
          is_computed_u(newind) = true;
          top = top + 1;
          list(top) = newind;
          is_row(top) = true;
        end
      end      
     end
    curr = curr + 1;
 end

 bestred = 0;
 found = false;
 for count = 1:m
    i = next_row;
    %for j = 1:n
    %  redcost = costm(i,j) - u(i) - v(j);
      [redcost,j]=min(costm(i,:) - u(i) - v'); % vectorized
      if redcost < bestred 
        bestred = redcost;
        indi = i;
        indj = j;
      end
    %end
    next_row=next_row+1;
    if (next_row > m)
      next_row = 1;
    end
    if (bestred < THRESHOLD)
      found=true;
      break;
    end
 end
  % indi and indj are
  % i and j index of new basic variable 
 %% DEBUG
 %found
 %disp([indi, indj])
 %costm-u-v'
 %assert(basis(indi,indj)==false,'Re-selected base element');
 %% 
 if (found)
   %% add_to_basis
  % disp('in addtobasis'); % DBG
  basis(indi, indj) = true;
  basis_byrow_top(indi) = basis_byrow_top(indi)+1;
  basis_byrow( indi, basis_byrow_top(indi)) = indj;

  basis_bycol_top(indj) = basis_bycol_top(indj)+1;
  basis_bycol( indj, basis_bycol_top(indj)) = indi;
  %% DEBUG
  %basis
 % basis_byrow
 % basis_byrow_top'
 % basis_bycol
 % basis_bycol_top'
 % assert(sum(basis(:))==m+n,'Augmented basis broken')
  %%
  %% find_circle
  %disp('in findcircle') %DBG
  circlea(:) = 0; %DBG
  circleb(:) = 0; %DBG
  candlist(:) = 0; %DBG
  circlea(1) = indi;
  circleb(1) = indj;
  curr = 1;           
  curr_fork = 0; %  -1;  
                      
  finished = false;

  do_rowscan = true;     
  next_branch = 1;
 % GFXDBG 
 % hold off;imagesc(basis);colormap(flipud(gray));hold on;
  while (~finished)
    lasti = circlea(curr);
    lastj = circleb(curr);
    candi = 0; candj = 0;
    %%DEBUG
    %[circlea(1:curr),...
    % circleb(1:curr)]'
  % GFXDBG
  %  plot(circleb(1:curr),circlea(1:curr),'-o'); pause();
    ncand = 0;
%    disp('in countcandidates');
    if (do_rowscan) 
     % count candidates and compile candlist 
      for j = 1:basis_byrow_top(lasti)
        candj = basis_byrow( lasti, j);      
        if (basis_bycol_top(candj) > 1 && candj ~= lastj) 
          ncand=ncand+1;
          candlist(ncand) = candj;
        end
      end
    else
   % COLSCAN 
      for i = 1:basis_bycol_top(lastj)
        candi = basis_bycol(lastj, i);
        if candi == indi && curr > 3
          finished = true;
          break;  % the for loop
        end
        if (basis_byrow_top(candi) > 1 && candi ~= lasti)
           ncand=ncand+1;
           candlist(ncand) = candi;
        end
      end
      if (finished) 
        break;  
      end     
    end
      %% DEBUG
      %    next_branch, ncand
      %    if(ncand>0), candlist(1:ncand)', end
      %if(next_branch==2 && curr_fork==0),
     %  fprintf('lst %d %d lst1 %d nb %d ncnd %d crr %d cfrk %d\n',...
     %     lasti, lastj, candlist(1), next_branch,ncand,curr,curr_fork);
      %end
   %% 
    if ncand == 0
    % dead end (leaf), try new branch from last fork
      curr = rem_curr(curr_fork); 
      do_rowscan = rem_do_rowscan(curr_fork);
      next_branch = rem_next_branch(curr_fork);
    elseif ncand == 1 
    % add unique choice to circle
      curr = curr+1;
      if(do_rowscan)
        circlea(curr) = lasti;
        circleb(curr) = candlist(1);
      else
        circlea(curr) = candlist(1);
        circleb(curr) = lastj;
      end
      do_rowscan = ~do_rowscan; % change direction
      next_branch = 1;
    else %  ncand >= 2
      if next_branch == 1      
       % push state to "remember stack"
        curr_fork = curr_fork + 1;
        rem_curr(curr_fork) = curr; % backup position
        rem_do_rowscan(curr_fork) = do_rowscan;
        rem_next_branch(curr_fork) = 2; 
          
        curr = curr + 1;  % add new candidate
        if(do_rowscan)
          circlea(curr) = lasti;
          circleb(curr) = candlist(1);
        else
          circlea(curr) = candlist(1);
          circleb(curr) = lastj;
        end 
          
        do_rowscan =~do_rowscan; 
      elseif (next_branch > ncand )
        % dead end (all branches starting at current fork)
        %   try next branch 
        curr_fork = curr_fork - 1; 
        curr = rem_curr(curr_fork);
        do_rowscan = rem_do_rowscan(curr_fork);
        next_branch = rem_next_branch(curr_fork);
      else
        rem_next_branch(curr_fork) = rem_next_branch(curr_fork) + 1;
        curr = curr + 1;
        if(do_rowscan)
          circlea(curr) = lasti;
          circleb(curr) = candlist(next_branch);
        else
          circlea(curr) = candlist(next_branch);
          circleb(curr) = lastj;
        end 
        do_rowscan = ~do_rowscan;
        next_branch = 1;
      end
    end % if 
  end % while
  circ_last = curr;
  %%DEBUG
  %[circlea(1:circ_last),...
  %circleb(1:circ_last)]'
  %assert(bitand(circ_last,1)==0,'circle does not contain even number');
  %%

  %% move_mass
  % minimize mass over "all minus signs"
  mass = assignment(circlea(2), circleb(2));
  argmin = 2;    
  for k = 4:2:circ_last
    newmass = assignment( circlea(k), circleb(k));
    if newmass < mass 
      mass = newmass;
      argmin = k;
    end
  end

  if mass > 0 
    for k = 1:2:circ_last 
      assignment(circlea(k), circleb(k)) = assignment(circlea(k), circleb(k)) + mass;
      assignment(circlea(k+1), circleb(k+1)) = assignment(circlea(k+1), circleb(k+1)) - mass;
    end
  %  disp(iter);
  %spy(assignment); drawnow% DBG
  end

  % index of basis vector which shall be removed */
  indi = circlea(argmin);
  indj = circleb(argmin);

  %%  remove_from_basis
  %% DEBUG
  %indi,indj
  %%
  %assert(basis(indi,indj),'Trying to remove non-basis element');
  basis( indi, indj) = false;

  if basis_byrow_top(indi)==1 
    basis_byrow_top(indi) = 0;
  else 
    for j = 1:basis_byrow_top(indi) 
      if basis_byrow( indi, j) == indj 
        basis_byrow( indi, j) = basis_byrow( indi, basis_byrow_top(indi));
        basis_byrow( indi, basis_byrow_top(indi))=0; % DEBUG
        basis_byrow_top(indi)= basis_byrow_top(indi)-1;
        break;
      end
    end
  end

  if basis_bycol_top(indj)==1
    basis_bycol_top(indj) = 0;
  else 
    for (i = 1:basis_bycol_top(indj)) 
      if (basis_bycol( indj, i) == indi) 
        basis_bycol( indj, i) = basis_bycol( indj, basis_bycol_top(indj));
        basis_bycol( indj, basis_bycol_top(indj))=0; % DEBUG
        basis_bycol_top(indj)= basis_bycol_top(indj)-1;
        break;
      end
    end
  end
  end
 end
 cost=assignment(:)'*costm(:);
end


function test
%%
C=[3 7 6 6 ; 1 6 8 8 ; 3 0 8 1; 0 7 9 9]
revsimplex(ones(4,1),ones(4,1),C)
%%
 b=[1,5,2,8,2]';
 a=[3,8,1,6]';
 C=[3 4 6 8 9; 2 2 4 5 5; 2 2 2 3 2; 3 3 2 4 2];
 [p,c]=revsimplex(a,b,C)
 [q,e]=revsimplex(b,a,C')
%%
C=[3 7 6 6 ; 1 6 8 8 ; 3 0 8 1; 0 7 9 9];
revsimplex(ones(8,1),2*ones(4,1),[C,C]')' % works
revsimplex(2*ones(4,1),ones(8,1),[C,C])  % fails
end
