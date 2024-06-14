function [delta,Si,acceptL,Seps]=deltamim(x,y,varargin)
%% DELTAMIM Estimation of Borgonovo's delta moment independent measure
%			 (kernel density estimator included)
% DELTA = DELTAMIM(X,Y) returns the moment independent measure for input
%         arguments X and output argument Y
% DELTA = DELTAMIM(X,Y,M) additionally defines the number of partitions
% DELTA = DELTAMIM(X,Y,[],'Gfx Title') additionally produces a figure.
% [DELTA,SI] = DELTAMIM(...) additionally returns the first order effects
% [DELTA,SI,AL] = DELTAMIM(...) gives acceptance levels AL for the KS test
%
%    M may also be a structure or key-value arguments containing the fields
%       PartitionSize [min(sqrt(n),48)]
%       QuadraturePoints [110]
%       KSLevel [0.95] 0: off, >0: cheap KS,  <0: full KS
%       ZeroCrossing { ['on'] 'off' 'positive' 'negative' 'test'}
%       ParameterNames { 'x_1' ... 'x_n' 'y'}
%       KDEstimator { ['cheap'] 'stats' 'diffusion'
%                    'hist' 'nearestneighbor' 'pilot' }
%       KDWidth { ['auto'] 'clone' or number }
%       KDShape {'normal' 'triangle' ['epanechnikov']
%                    'box','biweight'} }
%       Complement  { ['off'] 'on' }
%       OutputTrafo { ['off'] 'cdf' 'normal' function_handle }
%       ShowOpts    { ['off'] 'on' }
%		ShowSep 	{ ['on'] 'off' 'only' }

% DDD { [false] true } - experimental 3d plot option
% PartitionSplit { [] } - quantiles for splitting the partition

% Written by elmar.plischke@tu-clausthal.de
%         Plischke, Borgonovo, Smith:
%         "Global sensitivity measures from given data",
%         European Journal of Operational Research 226(3):536-550, 2013

[n,k]=size(x);

%% process options
% defaults
opts=struct('PartitionSize',min(ceil(n^(2/(7+tanh((1500-n)/500)))),48),...
	           'QuadraturePoints',110,'KSLevel',.95,...
               'ZeroCrossing','on','ParameterNames',[],...
               'KDEstimator','cheap','KDWidth','auto',...
               'Complement','off','SwitchXY','off',...
               'OutputTrafo','off',...
               'PlotCols',ceil(sqrt(k)),'PlotCmd',@plot,...
               'ShowOpts','off',...
			         'ShowSep','on',...
               'KDShape','epanechnikov','DDD',false,...
	             'PowerLoss',0,'PowerLossScale',0,...
               'PartitionSplit',[],...
               'GfxTitle','');
% argument handling
if(nargin>=3)
 M=varargin{1};
 if(~isempty(M))
  if(isnumeric(M))
  % partition size / gfx Title
   opts.PartitionSize=M;
   if(nargin==4)
    gfx=varargin{2}; else, gfx=[];
   end
  elseif isstruct(M)
  % option structure / gfx Title
   members=fieldnames(opts);
   for i=1:length(members)
    o=members{i};
    if isfield(M,o), opts.(o)=M.(o);end
   end
   if(nargin==4)
    gfx=varargin{2}; else, gfx=[];
   end
  else
  % arguments are key value pairs
   M=struct(varargin{:});
   members=fieldnames(opts);
   for i=1:length(members)
    o=members{i};
    if isfield(M,o), opts.(o)=M.(o);end
   end
   gfx=opts.GfxTitle;
  end
 end
end
%% Construct default parameter names
if isempty(opts.ParameterNames)
    pnames=cell(1,k+1);
    for i=1:k
        pnames{i}=['x_{' num2str(i) '}'];
    end
    pnames{k+1}='y';
    opts.ParameterNames=pnames;
end
if(length(opts.ParameterNames)==k), opts.ParameterNames{k+1}='y'; end
if(~isempty(opts.PartitionSplit))
 [pp,kk]=size(opts.PartitionSplit);
 if(kk~=1), error('Only column vector in PartitionSplit allowed'); end
 opts.PartitionSize=pp+1;
end
if strcmpi(opts.KDEstimator,'diffusion')
    % round to next power of 2
    opts.QuadraturePoints=2^nextpow2(opts.QuadraturePoints);
end
%% set kernel
global Kernel % dirty hack for passing Kernel through to the density estimator...
% these kernels have std. deviation 1
switch lower(opts.KDShape)
    case 'normal'
        Kernel=@(x)exp(-x.^2/2)/sqrt(2*pi);
    case 'triangle'
        Kernel=@(x)max(1-abs(x/sqrt(6)),0)/sqrt(6);
    % R2023b built-in kde() calls it parabolic
    case {'epanechnikov','parabolic'}
        Kernel=@(x)3/(4*sqrt(5))*max(1-(x.^2/5),0);
    case {'box','uniform'}
        Kernel=@(x)(1-abs(x/sqrt(3))>0)/(2*sqrt(3));
    case {'biweight','biquadratic'}
        Kernel=@(x)15/(16*sqrt(7))*max(1-(x.^2/7),0).^2;
    otherwise
        error('Unsupported kernel');
end
clear Kernel % avoid complaints about reinitialising global vars
%%
% compute critical value for KS statistics
kolmog=@(x,y)(x<4)*sqrt(2*pi)/x*sum(exp(-(1:2:35).^2*pi^2./(8*x^2)))...
            +(x>=4)*1.0-y;
if(opts.KSLevel~=0)
    % Complement mode only works with cheap version
    if strcmpi(opts.Complement,'on'), opts.KSLevel=abs(opts.KSLevel); end
    if opts.KSLevel<0
        opts.KSLevel=-opts.KSLevel;
        KStest=@(s,y,f)max(abs(cumtrapz(y,f))); % full version
    else
        KStest=@(s,y,f)s*0.5;                   % cheap version
    end
    for i=1:length(opts.KSLevel)
        KScrit(i)=fzero(@(x)kolmog(x,opts.KSLevel(i)),[0.001,2]);
    end
else
    KScrit=0;
end
if isscalar(KScrit), KScrit=KScrit*ones(1,k); end

%%
if(~strcmpi(opts.ShowOpts,'off'))
    opts,
    if(opts.KSLevel~=0) KStest, end
end
%%
 if(isa(opts.OutputTrafo,'function_handle'))
    y=opts.OutputTrafo(y);
    opts.OutputTrafo='off';
 end
[ys,indx]=sort(y);
switch(lower(opts.OutputTrafo))
    case { 'off','none' }
% compute the grid
        ymaxmin=ys([1,end]);rangey=diff(ymaxmin);
% the support of the estimates will be larger than the data range:
% the less samples the more the kernel variation
% let's account for mass outside the range...
        ysupp = ymaxmin+0.06*rangey*[-1;1];
        yy= linspace(ymaxmin(1)-.04*rangey,ymaxmin(2)+.04*rangey,...
                     opts.QuadraturePoints);
    case 'cdf'
    % transform to empirical cdf
	      ys=empcdf(ys);
	      y(indx)=ys;
        ysupp=[-0.06,1.06];
        yy=linspace(-0.04,1.04,opts.QuadraturePoints);
    case 'normal'
	      ncdf=-sqrt(2)*erfinv(1-2*empcdf(ys)); % norminv without stats
        y(indx)=ncdf;ys=ncdf;
    % same as for untransformed
        ymaxmin=ys([1,end]);rangey=diff(ymaxmin);
        ysupp = ymaxmin+0.06*rangey*[-1;1];
        yy= linspace(ymaxmin(1)-.04*rangey,ymaxmin(2)+.04*rangey,...
        opts.QuadraturePoints);
    case 'interpol'
        % interpolate quadrature points
        % ZeroCrossing has to be off
        % opts.ZeroCrossing='off';
        yy=interp1((2*(1:n)-1)/(2*n),ys,...
            linspace(0,1,opts.QuadraturePoints),'spline');
    case 'cdf-tight'
        % the next two test influence of support
        ys=empcdf(ys);
        y(indx)=ys;
        ysupp=[-0.02,1.02];
        yy=linspace(-0.01,1.01,opts.QuadraturePoints);
    case 'cdf-loose'
        ys=empcdf(ys);
        y(indx)=ys;
        ysupp=[-0.1,1.1];
        yy=linspace(-0.08,1.08,opts.QuadraturePoints);
    otherwise
        error('Unsupported output transformation');
end
ysupp2=yy([ 1 end]);
yy_=yy; % a backup copy

xr=zeros(n,k);
%% transform x into ranks and sort
consticator=zeros(1,k); % constant input indicator
for i=1:k
    [xx,indxx]=sort(x(:,i));
    consticator(i)=xx(end)==xx(1);
    xx(indxx)=empcdf(xx);
    xr(:,i)=xx(indx); % empirical cdf sorted for y
end

%% smoothed empirical pdf
if(ischar(opts.KDWidth))
    alfa=0;
else
    alfa=opts.KDWidth;
end
switch lower(opts.KDEstimator)
    case 'cheap'
            [f1,alfa]=kdest(y,yy,alfa);
    case 'pilot'
            [f1,alfa]=kdepilot(y,yy,alfa);
    case 'stats'
        if(alfa==0)
            f1=ksdensity(y,yy,'support',ysupp,'kernel',opts.KDShape);
        else
            f1=ksdensity(y,yy,'support',ysupp,'kernel',opts.KDShape,...
               'width',alfa);
        end
    case 'diffusion'
        if(alfa==0)
          [alfa,f1]=kde(ys,opts.QuadraturePoints,ysupp2(1),ysupp2(2));
        else
          [alfa,f1]=kde(ys,opts.QuadraturePoints,ysupp2(1),ysupp2(2),alfa);
        end
        f1=f1';f1(f1<0)=0;
    case 'hist'
        f0=histc(ys,yy);
        f1=f0./trapz(yy(:),f0);
    case 'nearestneighbor'
        if(alfa==0)
	        f1=densNN(ys,yy(:));
        else
 	        f1=densNN(ys,yy(:),alfa(1));
        end
        f1=f1./trapz(yy,f1);
    otherwise
        error('Unknown kernel density estimator');
end
%%
if(opts.PowerLoss~=0)
  gamma=opts.PowerLoss-1;
 % EfYgm1=mean(interp1(yy,f1,ys).^gamma);
  EfYgm1=trapz(yy,f1.^(gamma+1));
  if(opts.PowerLossScale==0)
   opts.PowerLossScale=2/gamma/(gamma+1);
  end
end
%%
M=opts.PartitionSize;
if(isempty(opts.PartitionSplit))
  segs=linspace(0,1,M+1);
else
  segs=[0,opts.PartitionSplit',1];
end
%% prepare gfx output
if(~isempty(gfx))
  cols=jet(M);
  clf
	layoutrows=ceil(k/opts.PlotCols);
	if(strcmpi(opts.ShowSep,'on')), layoutrows=layoutrows+1; end
	if(strcmpi(opts.ShowSep,'only')), layoutrows=0; end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=zeros(1,k); % moment independent delta measure
Sr=zeros(1,M);    % conditional separations
nr=zeros(M,1);    % number of realisations per subspample

%% % additionally compute first order effects via CR
if(nargout>=2)
 Si=zeros(1,k);
 Vyc=zeros(1,M);
 Ey=mean(y);
 Vy=var(y);
end
if(nargout>=3)
    Kr=zeros(1,M);      % Kolmogorov per partition
    acceptL=zeros(2,k); % Rejection level simple/full
end
%% parameters of interest
for i=1:k
 if(consticator(i)), continue; end % skip constant input
 if(~isempty(gfx))
  if layoutrows>0
    if(k>1)
	    subplot(layoutrows,opts.PlotCols,i);
    else
      subplot(layoutrows,1,i);
    end
    if(~opts.DDD)
      opts.PlotCmd(yy,f1,'k','LineWidth',2);
      ylabel('Density function');
    else
      plot3(yy,0*f1,f1,'k','LineWidth',2);
      zlabel('Density function');
    end
    xlabel([ opts.ParameterNames{k+1} ' given ' opts.ParameterNames{i} ]);
    title(gfx);
    hold on
  end
 end
%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for m=1:M
    if(~strcmpi(opts.Complement,'on'))
        yx = ys(xr(:,i)>=segs(m) & xr(:,i)<segs(m+1));
    else
        yx = ys(xr(:,i)<segs(m) | xr(:,i)>=segs(m+1)); % test complementary
    end
    nx = length(yx);
    nr(m)=nx;
    Sr(m)=0;
    if (nx>0) % non-empty segment
%%% contribution to the variance of the conditional expectation
      if(nargout>=2)
        if(~strcmpi(opts.Complement,'on'))
          Vyc(m)= nx*(mean(yx)-Ey)^2;
        else
          Vyc(m)= nx^2/(n-nx)*(mean(yx)-Ey)^2;
        end
      end
 %%%
      if(strcmpi(opts.KDWidth,'auto')), alfa=0; end
      switch lower(opts.KDEstimator)
      case 'cheap'
        f2=kdest(yx,yy,alfa);
      case 'pilot'
        f2=kdepilot(yx,yy,alfa);
      case 'stats'
        if(alfa~=0)
         f2=ksdensity(yx,yy,'support',ysupp,'kernel',opts.KDShape,...
                'width',alfa);
        else
         f2=ksdensity(yx,yy,'support',ysupp,'kernel',opts.KDShape);
        end
      case 'diffusion'
        if(alfa~=0)
         [~,f2]=kde(yx,opts.QuadraturePoints,ysupp2(1),ysupp2(2),alfa);
        else
         [~,f2]=kde(yx,opts.QuadraturePoints,ysupp2(1),ysupp2(2));
        end
        f2=f2';f2(f2<0)=0;
      case 'nearestneighbor'
        if(alfa==0)
         f2=densNN(yx,yy(:));
	      else
         f2=densNN(yx,yy(:),alfa(end));
        end
        f2=f2./trapz(yy,f2);
      case 'hist'
        f0=histc(yx,yy);
        f2=f0/trapz(yy(:),f0);
      end
      if(~isempty(gfx))
	   if(layoutrows>0)
        if(k>1)
          subplot(layoutrows,opts.PlotCols,i);
        else
          subplot(layoutrows,1,i);
        end
        if(~opts.DDD)
         if(~strcmpi(opts.Complement,'on'))
          opts.PlotCmd(yy,f2,'Color',cols(m,:));
         else
        % Complement mode
          opts.PlotCmd(yy,max(1/(n-nx)*(n*f1-nx*f2),0),'Color',cols(m,:));
         end
        else
           plot3(yy,m+0*f2,f2,'Color',cols(m,:));
        end
	   end
      end
    if(opts.PowerLoss~=0)
        Sr(m)=opts.PowerLossScale*(trapz(yy,f2.^(gamma+1))-EfYgm1);
    else
	  % compute differences
      fff=f1-f2;

      deltaoffset=0;
 % detect zero crossings (no impressive effect)
      switch lower(opts.ZeroCrossing)
      case 'on'
        ff=abs(fff); deltafactor=0.5;
        indx=find(fff(1:end-1).*fff(2:end)<0);

        if(~isempty(indx))
          yz=yy(indx)+ (yy(indx+1)-yy(indx)).*fff(indx)./(fff(indx)...
                                                      -fff(indx+1));
          for j=length(indx):-1:1 % inserting from the back
            l=indx(j);
            yy=[yy(1:l), yz(j), yy(l+1:end) ];
            ff=[ff(1:l), 0, ff(l+1:end) ];
            fff=[fff(1:l), 0, fff(l+1:end) ];
          end
        end
      case 'positive'
        ff=max(fff,0);deltafactor=1;
      case 'negative'
        ff=-min(fff,0);deltafactor=1;
      case 'off'
        ff=abs(fff); deltafactor=0.5;
      case 'test'
        ff=[ max(fff,0),-min(fff,0)];  deltafactor=0.5;
        yy=[ yy yy-yy(1)+yy(end)]; % separate positive and negative parts
      case 'min'
        ff=min(f1,f2);deltafactor=-1;
        deltaoffset=2;
      otherwise
        disp('Unknown ZeroCrossing method');
        ff=abs(fff); deltafactor=0.5;
      end

      if(strcmpi(opts.Complement,'on'))
        ff=nx/(n-nx)*ff;
        nr(m)=n-nx;
        nx=n-nx; % for the KS test below
      end

% integrate using trapzoidal rule
      S=deltaoffset+2*deltafactor*trapz(yy(:),ff(:));

% test for Kolmogorov-Smirnov (p=95%) (asymptotics for nx>30) (helps lots)
      if(KScrit(i) && KStest(S,yy,fff) < sqrt(1/n+1/nx)*KScrit(i))
        Sr(m)=0;
        if(nargout>=2), Vyc(m)=0; end % clear first order contribution
      else
        Sr(m)=S;
      end
% save KS for endogeneous threshold computation
      if(nargout>=3)
        Kr(m)=max(abs(cumtrapz(yy,fff)));
      end
      yy=yy_; % clear zero crossings
    end % ~PowerLoss
    end % non-empty segment
end

if(nargout>=3)
% predict rejection level
  thres1=max(.5*Sr./sqrt(1/n+1./nr') );
  thres2=max(Kr./sqrt(1/n+1./nr') );
  thres3=min(Kr./sqrt(1/n+1./nr') );
  if(thres1>0), acceptL(1,i)=kolmog(thres1,0); end
  if(thres2>0), acceptL(2,i)=kolmog(thres2,0); end
  if(thres3>0), acceptL(3,i)=kolmog(thres3,0); end
end

if(~isempty(gfx))
   if(layoutrows>0)
   % over-plot total output density
   if(~opts.DDD)
        opts.PlotCmd(yy,f1,'k','LineWidth',2);
        hold off
   else
       plot3(yy,0*f1+M+1,f1,'k','LineWidth',2);
       hold off
   end
   end
   if ~strcmpi(opts.ShowSep,'off')
    if(layoutrows>0)
     subplot(layoutrows,1,layoutrows);
	  else
	   subplot(1,1,1);
    end
    xcols=get(gca,'ColorOrder');
    plot( cumsum(nr)/n-1/(2*M),Sr,...
      ...   plot( cumsum(nr)/n,Sr,...
             'Color',xcols(1+mod(i-1,length(xcols)),:) );
			     hold on
  %% Colorbar testing
  %colormap(cols);
  %a=axis
  %alpha(image([a(1)+(a(2)-a(1))/M,a(2)-(a(2)-a(1))/M],...
  %    .2*[a(3)+(a(4)-a(3))/M,a(4)-(a(4)-a(3))/M],1:M),.2)
  %axis(a)
  end
 end
% Expectation over all partition segments
delta(i)=0.5*(Sr*nr)/n;
Seps(i,:)=Sr;
% first order sensitivity index (biased estimator)
if(nargout>=2)
	Si(i)=sum(Vyc)/Vy/(n-1);
end

end

if(~isempty(gfx))
 % label the separation graphs
   if ~strcmpi(opts.ShowSep,'off')
    if(layoutrows>0)
     subplot(layoutrows,1,layoutrows);
    end
    ylabel('S_r');
    %xlabel('Partition Intervals');
    xlabel('Empirical cdf of inputs');
    title('Separation of Conditional Densities');

    if(KScrit)
        q=opts.PartitionSize;
        plot([0,1],KScrit(i)*[1,1]*sqrt((q+1)/n),'k:')
        legend({opts.ParameterNames{1:k},'cut-off'});
    else
        legend(opts.ParameterNames{1:k});
    end
    hold off
 end
end

end

function [est,h]=kdest(y,z,h)
% BOWMAN Kernel density estimator by Bowman Azzalini
% EST=BOWMAN(Y,Z,H) returns an estimate of the density
%     at the evaluation points Z given data Y and smoothness H
%%
global Kernel
 n=length(y);
 if(nargin<3 || h==0)
    m=median(y);
    s=min(std(y),median(abs(m-y))/0.675);
    h=s/(((3*n)/4)^(1/5));
 end
 %k=length(z);
 %W=repmat(z,n,1)-repmat(y,1,k);
 %est=mean(Kernel(W/h))/h;
 est=mean(Kernel(bsxfun(@minus,z,y)/h))/h;
end
%% experimental kernel density
function [est,h]=kdepilot(y,z,h)
% KDEPILOT Two step kernel density estimator by Hengartner
%%
global Kernel
 n=length(y);
 if(nargin<3 || h==0)
    m=median(y);
    s=min(std(y),median(abs(m-y))/0.675);
    h=s/(((3*n)/4)^(1/5));
 end
 h0=1.5*h;
 % implicit bsxfun
 alfa=mean(Kernel((y'-y)/h0))/h0;
 est0=mean(Kernel((z-y)/h0))/h0;
 est=((1./alfa)*(Kernel((z-y)/h))/h/n).*est0;
end

function xr = empcdf(xs)
% EMPCDF Empirical cdf with ties
  n=size(xs,1);
  xr=(1:n)';						 % ranks
  tie_loc=[find(diff(xs)==0);n+2]; % stopper
  tie_next=diff(tie_loc);
  maxt=numel(tie_loc);
  i=1;
  while(i < maxt)
	run=tie_loc(i);len=1;
	while(tie_next(i)==1), i=i+1;len=len+1; end
	xr( run : run + len) = run+len/2;
	i=i+1;
  end
  xr=(xr-.5)/n;
end

%% alternative kernel density estimator
% using built-in dct/idct
function [bandwidth,density,xmesh]=kde(data,n,MIN,MAX,bandwidth_in)
% Reliable and extremely fast kernel density estimator for one-dimensional data;
%        Gaussian kernel is assumed and the bandwidth is chosen automatically;
%        Unlike many other implementations, this one is immune to problems
%        caused by multimodal densities with widely separated modes (see example). The
%        estimation does not deteriorate for multimodal densities, because we never assume
%        a parametric model for the data.
% INPUTS:
%     data    - a vector of data from which the density estimate is constructed;
%          n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12;
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where Range=max(data)-min(data);
% OUTPUTS:
%   bandwidth - the optimal bandwidth (Gaussian kernel assumed);
%     density - column vector of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%             - If no output is requested, then the code automatically plots a graph of
%               the density estimate.
%  Reference: (please cite in your work)
% Z. I. Botev, J. F. Grotowski and D. P. Kroese
%  "Kernel Density Estimation Via Diffusion"
%  Annals of Statistics,2010, to appear

%
%  Example:
%           data=[randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%              kde(data,2^14,min(data)-5,max(data)+5);


%  Notes:   If you have a more reliable and accurate one-dimensional kernel density
%           estimation software, please email me at botev@maths.uq.edu.au

%%
data=data(:); %make data a column vector
if nargin<2 % if n is not supplied switch to the default
    n=2^14;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;
%%
if nargin<4 %define the default  interval [MIN,MAX]
    minimum=min(data); maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/10; MAX=maximum+Range/10;
end
%% set up the grid over which the density estimate is computed;
R=MAX-MIN; dx=R/(n-1); xmesh=MIN+(0:dx:R); N=length(data);
%bin the data uniformly using the grid define above;
initial_data=histc(data,xmesh)/N;
%%
a= dct(initial_data); % discrete cosine transform of initial data
a=a*sqrt(2*n); a(1)=a(1)/sqrt(2); % account for different scaling

if(nargin<5)
% now compute the optimal bandwidth^2 using the referenced method
I=[1:n-1]'.^2; a2=(a(2:end)/2).^2;
%a2=2*N*(a(2:end)/2).^2
% use  fzero to solve the equation t=zeta*gamma^[5](t)
    t_star=fzero(@(t)fixed_point(t,N,I,a2),[0,.1]); % [0,.1]
else
    t_star=(bandwidth_in/R)^2;
end
% smooth the discrete cosine transform of initial data using t_star
a_t=a.*exp(-(0:n-1)'.^2*pi^2*t_star/2);
% now apply the inverse discrete cosine transform
if (nargout>1)||(nargout==0)
    a_t=a_t/sqrt(2*n);a_t(1)=a_t(1)*sqrt(2);
    density=idct(a_t)/n/R;
end
% take the rescaling of the data into account
bandwidth=sqrt(t_star)*R;
if nargout==0
    figure(1), plot(xmesh,density)
end
end
%################################################################
function  out=fixed_point(t,N,I,a2)
% this implements the function t-zeta*gamma^[l](t)
l=7;
f=2*pi^(2*l)*sum(I.^l.*a2.*exp(-I*pi^2*t));
for s=l-1:-1:2
    K0=prod(1:2:2*s-1)/sqrt(2*pi);  const=(1+(1/2)^(s+1/2))/3;
    time=(2*const*K0/N/f)^(2/(3+2*s));
    f=2*pi^(2*s)*sum(I.^s.*a2.*exp(-I*pi^2*time));
end
out=t-(2*N*sqrt(pi)*f)^(-2/5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=densNN(xx,xs,k)
% DENSNN kth Nearest Neighbor Density

% Reference: G. Biau, L. Devroye: Lectures on the Nearest Neighbor Method,
% Springer 2015
[n,d]=size(xx);
[~,dd]=size(xs);
if(nargin<3), k=ceil(sqrt(n)); end
if(d~=dd), error('densNN: Dimension mismatch'); end
if any(k>n), error('densNN: Out of samples');end
dist=sort(bsxfun(@plus,sum(xx.*xx,2),sum(xs.*xs,2)')-2*xx*xs');
Vd=pi^(d/2)/gamma(d/2+1);

f=bsxfun(@rdivide,k(:),n*Vd*dist(k,:).^(d/2));
end
