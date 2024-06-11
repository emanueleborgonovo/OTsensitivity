function [perm,total_cost]=hungarian(a)
% Hungarian Kuhn Munkres Algorithm
% P=HUNGARIAN(A)
% A n x m cost matrix  
 
% Andrej Lopatin O(n^2 m)
% http://e-maxx.ru/algo/assignment_hungary

  [n,m]=size(a);
	
  vector=@(n)zeros(n,1); % ,'int32');
  u=vector(n); v=vector(m+1); p=vector(m+1); way=vector(m);
  
  for i=1:n 
	j0 = m+1; % marks the root (as 0 is no valid index in matlab)
	p(j0) = i;

	minv = inf(m,1);
	used = false(m+1,1);
	 %do
    do_first=true;
    while do_first || p(j0) 
	  used(j0) = true;
	  i0 = p(j0);  delta = inf; j1=j0;
	  for j=1:m
			if (~used(j)) 
				cur = a(i0,j)-u(i0)-v(j);
				if (cur < minv(j))
					minv(j) = cur;  way(j) = j0;
				end
				if (minv(j) < delta)
					delta = minv(j);  j1 = j;
				end
			end
		end
		for (j=m+1:-1:1) % include the root 
			if (used(j))
				u(p(j)) = u(p(j)) + delta;
                v(j)    = v(j)    - delta;
			else
				minv(j) = minv(j)-delta;
			end
		end
		% u(i) = u(i)+delta; % j=m+1
		j0 = j1;
	    do_first=false;
   end
   %until (p(j0) == 0) % free column
   do_first=true;
   %do
   while do_first || j0~=m+1
		j1 = way(j0);
		p(j0) = p(j1);
		j0 = j1;
        do_first=false;
   end
	 %until(j0==m+1)
  end
  total_cost = 0;
  for j = 1:m 
	 total_cost = total_cost + a(p(j),j);
  end
  perm(p(1:m))=1:m;
end
