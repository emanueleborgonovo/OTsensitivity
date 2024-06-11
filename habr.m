function [B,cost]=habr(C)
% The Linear Assignment Problem, Jaroslav Habr (1961) approximate solution
% as of German Wikipedia for "Ungarische Methode"
n=size(C,1);
assert(size(C,2)==n,'Square matrix required');

D=C-mean(C,1)-mean(C,2)+mean(C(:)); % subtract row and column mean, add overall mean

% build a latin cube (single entries per row/column) using the smallest values of D
[~,ij]=sort(D(:));
ii=1+mod(ij-1,n);
jj=1+(ij-ii)/n;

B=zeros(n,1);
marked=false(n,2);
k=1;c=0;
while c<n
 i=ii(k);j=jj(k);
 if(~marked(i,1) && ~marked(j,2))
  B(i)=j;
  marked(i,1)=true;
  marked(j,2)=true;
  c=c+1;
 end
 k=k+1;
end
 
 cost=0;for i=1:n, cost=cost+C(i,B(i)); end

end

function testxx
%%
% Vater vernimmt Streit aus dem Kinderzimmer. Seine vier Kinder 
% Anna, Berta, Chiara und David zanken sich wieder einmal um die Spielsachen: 
% Eisenbahn, Kaufmannsladen, Puppe und den Zoo. Da es zu keiner friedlichen 
% Einigung kommt, schreitet Vater ein und befragt die Kinder nach der Rangordnung 
% ihrer Vorlieben bezüglich der vier Spielzeuge. Aus diesen Rangordnungen 
% bildet Vater eine 4x4-Matrix und versucht, durch geschickte Zuordnung der 
% Spielzeuge zu den Kindern die "Summe der Tränen" zu minimieren.
%   A  B  C  D 
C=[ 1, 1 ,1 ,2; ...
    3, 2, 4, 1; ...
	4, 4, 2, 4; ...
	2, 3, 3, 3];
[perm,cost]=habr(C)
%%
end
