function u=sobolpoints(n,d)
% SOBOLPOINTS Quasi Monte Carlo Low Discrepancy Series, Matlab stub.
 u = net(sobolpoints(d),n+1);
 u(1,:)=[]; % delete first realization
end