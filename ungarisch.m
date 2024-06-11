function [a,p]=ungarisch(C,use_kuhnmunkres)
%% Ungarisch
% Hungarian Assignment Problem 
% from Date / Nagi
% SOMETIMES DEADLOCKS on float

if(nargin>1)
[a,b]=classical_hungarian(C);
else
[a,b]=alternating_tree_hungarian(C);
end

p=0;
a=a';
for i=1:length(a), p=p+C(i,a(i)); end
end
%% Kuhn Munkres
% The classical hungarian algorithm, covering all zeros in the reducted 
% cost C-Dr-Dc by selecting n rows and columns
function [Ar,Ac]=classical_hungarian(C)

n=size(C,1);
assert(n==size(C,2),'square matrix required');

Ar=-ones(n,1);
Ac=-ones(n,1);
ST=zeros(n,1); % stack
Z=cell(n,1);   % adjacency list
% initial reduction
Dr=min(C,[],2);
Dc=min(C-Dr,[],1);

Vr=false(n,1); % rows covered
Vc=false(n,1); % columns covered
Pr=-ones(n,1); % predecessor row
Pc=-ones(n,1); % predecessor column

top=0; % empty stack
    
while ~all(Vr) % optimality check
%% debugging
%    match_count
% %    if(top>0), ST(1:top), end
% %    disp('RELATIVE COSTS')
% %    round(C-Dr-Dc)
%    disp('rows'),disp(find(Vr)') % , Dr'  
%    disp('cols'),disp(find(Vc)') % , Dc
%    skip_update
%%
    for i=1:n
        if ~Vr(i)
            % ST.push(i)
            top=top+1;ST(top)=i;
        end
        % initialize adjacency list for row i
        Z{i}=[];
        for j=1:n
             if(abs(C(i,j)-Dr(i)-Dc(j))<10*max(eps,eps(C(i,j)))) % float compare
                 Z{i}= [ Z{i}, j ]; % Z[i].push(j) 
             end
        end
    end
    skip_update=false;
    while(top>0)
        i=ST(top);top=top-1;  %ST.top ST.pop
        while ~isempty(Z{i})
            j=Z{i}(1);        % Z[i].front
            Z{i}=Z{i}(2:end); % Z[i].pop
            inew=Ac(j);
            if(inew==i), continue; end % on next j
            if ~Vc(j) % if column is uncovered
                Pc(j)=i; % update predecessor index
                if inew==-1 % unassigned column
                    % augment the current assignments by 1
                    ccur=j;
                    while ccur~=-1 % until current row has no predecessor
                        rcur=Pc(ccur);
                        Ar(rcur)=ccur;
                        Ac(ccur)=rcur;
                        ccur=Pr(rcur); % update current column index
                    end
                    % end augment
                    skip_update=true; 
					top=0; % empty stack
					break
                else
                    top=top+1;ST(top)=inew;
                    Pr(inew)=j;     % update predecessor index
                    Vr(inew)=false; % uncover the row
                    Vc(j)=true;     % cover the column
                end
            end
        end
        % if(skip_update), break; end
    end
    if(~skip_update)
    % update the dual variables
     thta=min(min(C(~Vr,~Vc)-Dr(~Vr)-Dc(~Vc)));
     Dr=Dr+thta*(1/2-Vr);
     Dc=Dc+thta*(1/2-Vc');
    else
    % reset
     Vr=Ar~=-1;
     Vc=false(n,1);     
     Pr=-ones(n,1);
     Pc=-ones(n,1);
%     top=0; % empty stack
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%% alternating tree variant proposed by Lawler
function [Ar,Ac]=alternating_tree_hungarian(C)
n=size(C,1);ns=1:n;
assert(n==size(C,2),'square matrix required');

Ar=-ones(n,1);
Ac=-ones(n,1);
ST=zeros(n,1); % stack

% initial reduction
Dr=min(C,[],2);
Dc=min(C-Dr,[],1);

%
precision= 10*max(eps(max(C,[],2)),eps);

% reset 
Vr=false(n,1); 
Vc=false(n,1);     
Pr=-ones(n,1);
Pc=-ones(n,1);
slack=inf(1,n); % a row
top=0; % empty stack
     
while ~all(Vr) % optimality check 
    skip_update=false; 
%
    % Augmenting path search
    ii=find(~Vr);topadd=length(ii);
    ST(top+(1:topadd))=ii; top=top+topadd; % bulk initialize stack
    while(top>0)
        i=ST(top); % ST.top 
        top=top-1; % ST.pop
        
        credr = max(C(i,:)-Dr(i)-Dc,0); % DEBUG enforce nonnegative
        ii=find(slack>credr);
        Pc(ii)=i;
        slack(ii)=credr(ii);
        
        %jj=credr<=10*max(eps(C(i,:)),eps); % ==0
        jj=find(credr<=precision(i)); % ==0
         slack(jj)=0;          % DEBUG enforce 0
         Dc(jj)=C(i,jj)-Dr(i); % DEBUG enforce 0
        for j_iter=1:length(jj)
            j=ns(jj(j_iter));
        %for j=ns(jj)
                inew=Ac(j);
                %if(inew==i),disp('continue');continue; end % on next j (missing in article?)
                if ~Vc(j) % if column is uncovered
                    if inew==-1 % unassigned column
                    % augment the current assignments by 1
                        ccur=j;
                        while ccur~=-1 % until current row has no predecessor
                            rcur=Pc(ccur);
                            Ar(rcur)=ccur;
                            Ac(ccur)=rcur;
                            ccur=Pr(rcur); %update current column index
                        end
                    % end augment
                        skip_update=true; break % out of for
                    else
                        top=top+1;ST(top)=inew;
                        Pr(inew)=j;     % update predecessor index
                        Vr(inew)=false; % uncover the row
                        Vc(j)=true;     % cover the column
                    end
                end
            %end
        end
        if(skip_update), break; end
    end
    if(~skip_update)
    % update the dual variables
     ss=slack>0; % >10*eps; % 0;
   %  if(isempty(ss)), disp('no tight slack?'); end
     thta=min(slack(ss));
   %  if(~isfinite(thta)), disp('theta infinite?'); end
     Dr=Dr+thta*(1/2-Vr); %(~Vr-Vr)/2;
     Dc=Dc+thta*(1/2-Vc)';
     slack(ss)=slack(ss)-thta;
     jj=ss & (slack<=16*max(eps,eps(thta))); %==0;
     % jj(j)=true; % DEBUG 
     topadd=sum(jj); % should never be empty
     slack(jj)=0; % DEBUG
     ST(top+(1:topadd))=Pc(jj); top=top+topadd; % bulk initialize stack
     skip_update=true;
    else 
    % % reset
     Vr=Ar~=-1;
     Vc=false(n,1);     
     Pr=-ones(n,1);
     Pc=-ones(n,1);
     slack=inf(1,n);
     top=0; % empty stack
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
function test
%%
C=[3 7 6 6 ; 1 6 8 8 ; 3 0 8 1; 0 7 9 9]
ungarisch(C)
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%
function test2
%%
ishigami
n=8192;
x=trafo(sobolpoints(n,k));
y=model(x);

[~,ix]=sort(x);

nextplot=@(i)subplot(1,2,i);
cnt=1;
C=round(-2*(y*y')+(sum(y.^2,2)+sum(y.^2,2)')); % integer Euclidean squared

type='hungarian';
M=44;
ms=round(linspace(0,n,M+1));
W=zeros(M,k);
tic
for m=1:M
   
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices
     
    for i=1:k
    % conditional output
    
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     [a,p]=ungarisch(C(jj,ii)); % linear assignment problem
     W(m,i)=p/nn;
    end
end
disp('Hungarian algorithm (downscale)');t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('hungarian (downscale)');
%end
%if bitget(selection,6) % currently
type='hungarianupscale';
W=zeros(M,k);
tic
for m=1:M
    m
     mc=ms(m+1)-ms(m);
     jj=round(linspace(.5,mc+.499,n)); % upscaling indices
     
    for i=1:k
    % conditional output
     i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     % upscale and compute
     [a,b]=ungarisch(C(:,ii(jj))); % linear assignment problem
     W(m,i)=b/n;
    end
end
disp('Hungarian algorithm (upscale)');t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('hungarian (upscale)');

%%
end

%
function failed
%%
for r=1:5
 r
 x=rand(40,3);
 C=-2*x*x'+sum(x.^2,2)+sum(x.^2,2)';
 D=C(repmat(1:10,[1,4]),:);

 tic,[p,c]=transsimp2(ones(40,1),ones(10,1)*4,C(1:10,:)');c,toc
 tic,[p,c]=transsimp2(ones(10,1)*4,ones(40,1),C(1:10,:));c,toc
 tic,[p,c]=ungarisch(D,1);c,toc
 
 tic,[p,c]=ungarisch(D);c,toc
 tic,[p,c]=ungarisch(D');c,toc
 
 tic;[p,c]=lapjv(D);c,toc
end 
%%
end