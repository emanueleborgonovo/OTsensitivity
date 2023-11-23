function [results]=bigwassersteintest(x,y,M,selection)
%% big Wasserstein Test
if nargin<4 || isempty(selection)
 selection = hex2dec('FFFF'); % all bits set
end
results=struct();

modenames={'simplex','sinkhorn','swap','hungarian','habr','lawler','bures',...
    'gradient','divergence','sinkscale','mack','revsimplex','otpython','sinkmem'};
needcost=[1 1 0 1 1 1 0 0 1 1 1 1 1 0];
MAXBIT=length(modenames); % row*col>= MAXBIT fopr gfx
if(iscellstr(selection))
 modes=0;
 for str = selection
     bit=find(strcmp(str,modenames));
     if ~isempty(bit)
     modes=bitset(modes,bit);
     else
         fprintf('wasserstein: unknown mode %s\n',str{1});
     end
 end
 selection=modes;
end
l=bitcount(selection);
switch(l)
  case 1,
    nextplot=@(cnt)drawnow();
  case 2,
    nextplot = @(cnt)subplot(1,2,cnt);
  case 3,
    nextplot = @(cnt)subplot(1,3,cnt);
  case 4,
    nextplot = @(cnt)subplot(2,2,cnt);
  case {5,6},
    nextplot = @(cnt)subplot(2,3,cnt);
  case {7,8,9},
    nextplot = @(cnt)subplot(3,3,cnt);
  otherwise,
    nextplot = @(cnt)subplot(4,4,cnt);
end
cnt=1;
[n,k]=size(x);
if(any(needcost(logical(bitget(selection,1:MAXBIT)))))
 disp('cost matrix construction');
 tic
 C=-2*(y*y')+(sum(y.^2,2)+sum(y.^2,2)'); % Euclidean squared
 toc
 else
 disp('method(s) use(s) builtin costs');
 end
 ms=round(linspace(0,n,M+1));
[~,ix]=sort(x);
%%
if bitget(selection,1)
type=modenames{1};
disp('transport simplex');
W=zeros(M,k);
tic
for m=1:M
    m
     for i=1:k
     %    i
         ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
        % [~,S]=transsimp2(ones(n,1)/n,ones(nn,1)/nn,C(:,ii));
        [~,S]=transsimp2(ones(nn,1)/nn,ones(n,1)/n,C(ii,:));
         W(m,i)=S;
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Simplex');drawnow
end
%%
%%
if bitget(selection,2)
type=modenames{2};
disp('Sinkhorn');
W=zeros(M,k);
tic
for m=1:M
    m
     for i=1:k
      %   i
         ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
        [S,R]=sinkfastCost(-0.001,C(ii,:)); % 400 for environ, 2 for normal
       % [S,R]=sinkfastCost(100,C(ii,:));
         W(m,i)=S;
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;

nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Sinkhorn');drawnow
end
if bitget(selection,3)
type=modenames{3};
disp('Puccetti swap (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
    % i
     yc=y(ix(ms(m)+1:ms(m+1),i),:);
     jj=randperm(n,size(yc,1));
     % downscale and compute
     ii=wasserstein_swap2(-y(jj,:),yc,1);
     W(m,i)=mean(sum((yc(ii,:)-y(jj,:)).^2,2)); % nn*n/nn
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Swap*');drawnow
end
%{
if bitget(selection,4)
type='swapupscale';
disp('swap (upscale)');
W=zeros(M,k);
tic
for m=1:M
    m
     mc=ms(m+1)-ms(m);
     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
     %i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     % upscale and compute
     kk=wasserstein_swap2(-y(ii(jj),:),y,1);
     W(m,i)=mean(sum((y(kk,:)-y(ii(jj),:)).^2,2));
     %[a,b]=lapjv(C(:,ii(jj))); % linear assignment problem
     %W(m,i)=b/n;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('swap (upscale)');drawnow
end
%}
if bitget(selection,4)
type=modenames{4};
disp('Hungarian algorithm (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
     %i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     [a,b]=hungarian(C(jj,ii)); % linear assignment problem
     W(m,i)=b/nn;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Hungarian*');drawnow
end
if bitget(selection,5)
%type='jonkervolgenant';
%disp('Hungarian/Jonker Volgenant algorithm (downscale)');
type=modenames{5};
disp('Habr approximation (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
     %i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     %[a,b]=lapjv(C(jj,ii)); % linear assignment problem
	 [a,b]=habr(C(jj,ii)); % linear assignment problem
     W(m,i)=b/nn;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
%plot(1:M,W);title('Jonker/Volgenant*');drawnow
plot(1:M,W);title('Habr*');drawnow
end
%{
if bitget(selection,6)
type='hungarianupscale';
disp('Hungarian (upscale)');
W=zeros(M,k);
tic
for m=1:M
    m
     mc=ms(m+1)-ms(m);
     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
    % i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     % upscale and compute
     [a,b]=lapjv(C(:,ii(jj))); % linear assignment problem
     W(m,i)=b/n;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('hungarian (upscale)');drawnow
end
%}
if bitget(selection,6)
type=modenames{6};
disp('Hungarian/Lawler algorithm (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
     %i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     [a,b]=ungarisch(C(jj,ii));  % Lawler tree-based approach
     W(m,i)=b/nn;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Lawler*');drawnow
end
if bitget(selection,7)
    type=modenames{7};
    disp('Bures');
    tic;
    [~,W]=bwsi(x,y,M);
    t=toc
    results.(['W' type])=W;
    results.(['T' type])=t;
    nextplot(cnt);cnt=cnt+1;
    plot(1:M,W);title('Bures');drawnow
end
if bitget(selection,8)
    type=modenames{8};
    disp('Gradient');
    tic;
    W=wassgrad2(x,C,M,-0.001); % relative precision
 t=toc
 results.(['W' type])=W;
 results.(['T' type])=t;

 nextplot(cnt);cnt=cnt+1;
 plot(1:M,W);title('Gradient');drawnow
end
if bitget(selection,9)
type=modenames{9};
disp('Sinkhorn Divergence');
W=zeros(M,k);
tic
FF=[]; % reuse W(y,y) if it becomes available
for m=1:M
    m
     for i=1:k
      %   i
         ii=ix(ms(m)+1:ms(m+1),i);
         [w,opts]=sinkhorn4([],.05,y,y(ii,:),[],[],FF);
         W(m,i)=w;FF=opts.PotentialFF;
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;

nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Divergence');drawnow
end
%{
if bitget(selection,10)
type='simplexcomplement';
disp('transport simplex complement');
W=zeros(M,k);
tic
for m=1:M
    m
     for i=1:k
       %  i
         ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
        % [~,S]=transsimp2(ones(n,1)/n,ones(nn,1)/nn,C(:,ii));
        D=C(ii,:);D(:,ii)=[];
        [~,S]=transsimp2(ones(nn,1)*(1/nn-1/n),ones(n-nn,1)/n,D);
         W(m,i)=S;
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Simplex complement');drawnow
end
%}
if bitget(selection,10)
type=modenames{10};
disp('Sinkhorn (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
    % i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     a=sinkfastCost(-0.001,C(ii,jj));
     W(m,i)=a/nn;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Sinkhorn*');drawnow
end
if bitget(selection,11)
type=modenames{11};
disp('Bradford / Mack (downscale)');
W=zeros(M,k);
tic
for m=1:M
    m
%     mc=ms(m+1)-ms(m);
%     jj=round(linspace(.5,mc+.499,n)); % upscaling indices

    for i=1:k
    % conditional output
    % i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     [a,b]=lapmack(C(ii,jj)); % Bradford / Mack min cost flow
     W(m,i)=b/nn;
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Mack*');drawnow
end
%{
if bitget(selection,12)
type='greenkhorn';
W=zeros(M,k);
tic
for m=1:M
    for i=1:k
    % conditional output
    % i
     ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
     jj=randperm(n,nn);
     % downscale and compute
     a=greenkhorn(exp(-C(ii,jj)/2),ones(nn,1),ones(1,nn),1000,false,...
         [],[]); % Greedy Sinkhorn
     W(m,i)=sum(a.*C(ii,jj),'all')/nn;
    end
end
disp('Greenkhorn (downscale)');t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('greekhorn (downscale)');drawnow
end
%}
if bitget(selection,12) % 0
type=modenames{12};
disp('Revised Simplex');
W=zeros(M,k);
tic
for m=1:M
    m
     for i=1:k
     %    i
        ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
        % [~,S]=transsimp2(ones(n,1)/n,ones(nn,1)/nn,C(:,ii));
        [~,S]=revsimplex(ones(nn,1)*n,ones(n,1)*nn,C(ii,:));
        %[~,S]=revsimplex(ones(n,1)*nn,ones(nn,1)*n,C(:,ii)); % very slow
        W(m,i)=S/(n*nn);
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Revised simplex');drawnow
end
if bitget(selection,13) % 0
type=modenames{13};
disp('EMD (python)');
if ispc
%if isempty(pyenv().Version)
% % setenv('CONDA_DLL_SEARCH_MODIFICATION_ENABLE',1) % ?? maybe powershell instead of cmd issues ??
% pyenv('Version','C:\\ProgramData\\Anaconda3\\pythonw.exe') % instead of plain python

% matlab central
 py_root_to_use=fileparts('C:\\ProgramData\\Anaconda3\\_conda.exe');
 assert(exist(py_root_to_use,'dir')==7,'Wrong Python setup. Please modify.')
 env=strsplit(getenv('PATH'),';');
 add_to_path={
     fullfile(py_root_to_use, 'Library','mingw-w64','bin')
     fullfile(py_root_to_use, 'Library','usr','bin')
     fullfile(py_root_to_use, 'Library','bin')
     fullfile(py_root_to_use, 'Scripts')
     };
 env = unique([add_to_path(:); env(:)],'stable');
 setenv('PATH',strjoin(env,';'));
% OpenMP issues with mkl:
% As an unsafe, unsupported, undocumented workaround you can
% set the environment variable
 setenv('KMP_DUPLICATE_LIB_OK','TRUE');
% to allow the program to continue to execute
end
otmod=py.importlib.import_module('ot');
W=zeros(M,k);
tic
for m=1:M
    m
     for i=1:k
     %    i
        ii=ix(ms(m)+1:ms(m+1),i);nn=size(ii,1);
		S=double(otmod.emd2(py.numpy.array(ones(nn,1)*n),...
		                    py.numpy.array(ones(n,1)*nn),...
							py.numpy.array(C(ii,:))));
        W(m,i)=S/(n*nn);
     end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Python OT');drawnow
end
if bitget(selection,14) % 0
 type=modenames{14};
 disp('Sinkhorn (small memory footprint)');
W=zeros(M,k);
tic
for m=1:M
    m
    for i=1:k
    % conditional output
     yc=y(ix(ms(m)+1:ms(m+1),i),:);
     W(m,i)=sinkfast(-.001,y,yc);
     %W(m,i)=sinkfast(.05,y,yc);
    end
end
t=toc
results.(['W' type])=W;
results.(['T' type])=t;
nextplot(cnt);cnt=cnt+1;
plot(1:M,W);title('Sinkhorn (small memory)');drawnow
end
end % function
function testbig
%%
n=1024;
k=4;
u=net(scramble(sobolset(k),'MatousekAffineOwen'),n);
xx=u.*[13-7 0.12-0.02 3-0.01 30.295-30.01] + [7 0.02 0.01 30.01];
names={'Ma','Di','Lo','\tau'};
s = [0.5, 1, 1.5, 2, 2.5];
t = [0.3:0.3:60];
yy=zeros(n,length(s)*length(t));
tic
for i=1:n
    y= environ(xx(i,:));
 yy(i,:)=y(:)';
end
toc
res=bigwassersteintest(xx,yy,10,63);
 %%
Wsimplex=mean(sqrt(res.Wsimplex))
WSink=mean(sqrt(res.Wsinkhorn))
Wswap=mean(sqrt(res.Wswap))
WHung=mean(sqrt(res.Whungarian))
Wbures=mean(sqrt(res.Wbures))

Wsimplex/max(Wsimplex)
WSink/max(WSink)
Wswap/max(Wswap)
WHung/max(WHung)
Wbures/max(Wbures)
%%
%% and now for Gaussian
W2analyt=[6.47 6.52 2.86]; % ??
%%
r12=0.5;
r13=0.5;
r23=0.5;
muX=[1 1 1]
SigX=[1 r12 r13; r12 1 r23; r13 r23 1];
b=[0 0]';
A=[4 -2 1; 2 5 -1];
X=mvnrnd(muX,SigX,1000);
Y=(A*X'+b)';
lst=bigwassersteintest(X,Y,10,63);
%%
%results: first row analytical')
R=[W2analyt;sqrt(mean(lst.Wsimplex));sqrt(mean(lst.Wsinkhorn));...
    sqrt(mean(lst.Wswap));sqrt(mean(lst.Whungarian));sqrt(mean(lst.Whungarianupscale));sqrt(mean(lst.Wbures))]
%subplot(2,3,6)
%plot(1:5,R(2:end,:),'*',[1;5],[1;1]*R(1,:),'k:')
end

%% moretests
function washishi
%%
ishigami
for n=[1024 2048 4096]
 x=trafo(sobolpoints(n,k));
 y=model(x);
 for M=[8 16 32 64]
%simplex methods are [ 1 5 10 11 12]
% if M==8
% else
%modes=[3:6 11]; % all linear assignment problems / random selection
  modes =[2 3 4 5 7 10 11 12 13]; % jv, lawler too slow, sinkhorn, bures, simplex for comparison
% end
%dbstop if error % put this in the console
 diary on
 fprintf('%d/%d\n',n,M);
 diary off
 W=bigwassersteintest(x,y,M,sum(2.^(modes-1)));
 diary on
 W
 diary off
 end
end
%%
end
