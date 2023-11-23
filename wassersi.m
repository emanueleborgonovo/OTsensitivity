function [W,d,nn]=wassersi(x,y,M,gfx)
% WASSERSI Sensitivity indices using 1D Wasserstein distances
% W=WASSERSI(X,Y) returns a structure containing
% Wasserstein sensitivities: W1 , W2 , Winf
% squared Wasserstein sensitivity: W22
%
% W1 is bounded by 2MEAN(ABS(Y))
% W22 is bounded by 2VAR(Y))

% written by elmar.plischke@tu-clausthal.de
[n,k]=size(x);
[xs,ix]=sort(x);
ys=sort(y);
if(nargin<3) || isempty(M), M=10; end
if(nargin<4) || isempty(gfx), gfx=''; end

if(~isempty(gfx))
clf
L=sqrt(k);
if(ceil(L)*floor(L)>=k), myround=@floor; else myround=@ceil; end
cols=jet(M);
lso=get(gcf,'DefaultAxesLineStyleOrder');
set(gcf,'DefaultAxesLineStyleOrder',{'-'}); % just lines
end

Ms=round(linspace(1,n+1,M+1));
d=zeros(M,k);
e=zeros(M,k);
f=zeros(M,k);
g=zeros(M,k);
for m=1:M
  ii=Ms(m):Ms(m+1)-1;
  N=Ms(m+1)-Ms(m);
  nn(m)=N;
  jj=round(linspace(.5,N+.499,n));
  for i=1:k
     yc=sort(y(ix(ii,i)));
%      plot(xs(ii,i),yc,'ro',x(:,i),y,'b.');
% scale yc up to size n
     ycc=yc(jj);
     d(m,i)=mean((ys-ycc).^2);
	 e(m,i)=mean(abs(ys-ycc));
	 f(m,i)=max(abs(ys-ycc));
	 g(m,i)=(mean((ys-ycc).^8))^(1/8);
     if(~isempty(gfx))
     subplot(myround(L),ceil(L),i);plot(ys,ycc,'Color',cols(m,:),...
         'LineWidth',2); hold on
     xlabel('Output');ylabel('Conditional Outputs');
     title(['x_{' num2str(i) '}']);
     end
  end
end
if(gfx)
    set(gcf,'DefaultAxesLineStyleOrder',lso); %reset
end
W=struct('W2',nn*sqrt(d)/n,'W22',nn*d/n,'W1',nn*e/n,'W8',nn*g/n,'Winf',nn*f/n);
end

function testwassersi
%%
%for conditional normal
%    W2=(mY-mY|X)^2+(sigmaY-sigmaY|X)^2
%%
n=12300;
k=3;
model=@(x)x*[4;-2;1];
trafo=@(u)norminv(u,0,1);
x=trafo(rand(n,k));
y=model(x);
wassersi(x,y,10)
wassersi(x,y,20)
wassersi(x,y,50)
%% bring out the combine harvester
opt=struct('GfxTitle','blubb');
l=mydoubleloop(k,[64 64],model,trafo,[],opt);
l([7,6,8],:)
%% (semi) analytical solution
z=x(:,1);
Wana=[mean(sqrt(16*z.^2+(sqrt(21)-sqrt(5)).^2)), ...
    mean(sqrt(4*z.^2+(sqrt(21)-sqrt(17)).^2)), ...
    mean(sqrt(z.^2+(sqrt(21)-sqrt(20)).^2))]
%%
end