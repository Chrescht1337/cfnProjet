
function [L, U, P, Q] = LUsteps(A,pivot,simulaterounding);

% usage: [L, U, P, Q] = LUsteps(A,pivot,simulaterounding);
%    or: [L, U, P, Q] = LUsteps(A,pivot);
%    or: [L, U, P, Q] = LUsteps(A,simulaterounding);
%    or: [L, U, P, Q] = LUsteps(A);
%  input
%    A                   matrix
%    pivot               integer 0 (default),1,2 
%                        (no pivot, partial pivot, full pivot)
%    simulaterounding    real number (default 0) of simulated machine precision
%                        (0 for no rounding)
%  output
%    L, U                lower and upper triangular matrix
%    P, Q                row and column permutation matrices, so that
%                        so that P*A*Q = L*U

if nargin < 3, simulaterounding = 0; end
if nargin < 2, pivot = 0; end
if ~ismember(pivot,(0:2)), simulaterounding = pivot; pivot = 0; end
eps = simulaterounding;
[n m] = size(A);
P = eye(n); Q = eye(m);
nm = min(n,m);
L = [eye(nm); zeros(n-nm,nm)]; U = A;
% Gauss with pivots
for r = 1:nm,
   k = r; l = r;
   if pivot==1,
      [u k] = max(abs(U(r:n,r))); k = r-1+k;
   elseif pivot==2,
      [u k] = max(abs(U(r:n,r:n))); [u l] = max(u); k = k(l);
      k = r-1+k; l = r-1+l;
   end
   if k~=r, 
      % Ur = U(r,r:m); U(r,r:m) = U(k,r:m); U(k,r:m) = Ur; 
      % Lr = L(r,1:r-1); L(r,1:r-1) = L(k,1:r-1); L(k,1:r-1) = Lr;
      U([r k],1:n) = [0 1;1 0]*U([r k],1:n);
      P([r k],1:n) = [0 1;1 0]*P([r k],1:n);
      L([r k],1:r-1) = [0 1;1 0]*L([r k],1:r-1);
   end
   if l~=r,
      % Ur = U(1:n,r); U(1:n,r) = U(1:n,l); U(1:n,l) = Ur; 
      U(1:n,[r l]) = U(1:n,[r l])*[0 1;1 0];
      Q(1:m,[r l]) = Q(1:m,[r l])*[0 1;1 0];
   end
   rho = (rand(n-r,1)-0.5)*eps;
   L(r+1:n,r) = U(r+1:n,r)/U(r,r).*(1+rho);
   rho1 = (rand(n-r,m-r)-0.5)*eps;
   rho2 = (rand(n-r,m-r)-0.5)*eps;
   U(r+1:n,r+1:m) = (U(r+1:n,r+1:m)-(L(r+1:n,r)*U(r,r+1:m)).*(1+rho2))...
                     .*(1+rho2);
   U(r+1:n,r) = 0;
end
