prompt='Entrez la prÃ©cision machine';
eps=input(prompt);
matricule = 394691;
rand('state',matricule);
if exist('condA') ~= 1, condA = 1.e3; end
if exist('sizeA') ~= 1, sizeA = 100; end
n = sizeA;

% create matrices of size n x n with given condition number
% First make random matrix
X = rand(n,n);
% then find QR-factorisation to obtain an orthogonal matrix U
[U R] = qr(X);
% repeat the procedure to obtain an orthogonal matrix V
X = rand(n,n);
[V R] = qr(X);
% construct a diagonal matrix with condition number as given
Sigma = diag(1+(n-1:-1:0)/(n-1)*(condA-1));
A = U*Sigma*V';
% check that condition number is as given
condAtest = cond(A) 

x = floor(10*rand(n,1));
b = A*x;

eps2 = 1.e-8; % simulated machine precision (eps/2)
nsimul = 100; % number of simulation runs


% (1) condition
rand('state',matricule)
kappa = zeros(nsimul,1);
for simul = 1:nsimul,
   rho = (rand(n,1)-0.5)*eps2;

   bhat = b.*(1+rho);
   relerrbhat = norm((bhat-b)/b);%mean(abs((bhat-b)./b));%%
   xhat = A\bhat;
   relerrxhat = norm((xhat-x)/x);%mean(abs((xhat-x)./x));%%
   kappa(simul) = relerrxhat/relerrbhat;
end
plot(kappa)
kappacond = mean(kappa)/condA

% (2) stability
rand('state',matricule);
relerrxhat0 = zeros(nsimul,1);
relerrxhat1 = zeros(nsimul,1);
relerrxhat2 = zeros(nsimul,1);
for simul = 1:nsimul,
   [Lhat0 Uhat0] = LUsteps(A,eps2);
   y0=Lhat0\b;
   xhat0=Uhat0\y0;
   relerrxhat0(simul) = norm(abs((xhat0-x)/x));
   
   [Lhat1 Uhat1 P] = LUsteps(A,1,eps2);
   y1=Lhat1\(P*b);
   xhat1=Uhat1\y1;
   relerrxhat1(simul) = norm(abs((xhat1-x)/x));
   
   [Lhat2 Uhat2 P Q] = LUsteps(A,2,eps2);
   y2=Lhat2\(P*b);
   z=Uhat2\y2;
   xhat2=Q*z;
   relerrxhat2(simul) = norm(abs((xhat2-x)/x));
end
figure
plot(relerrxhat0,'b')
hold on
plot(relerrxhat1,'m')
plot(relerrxhat2,'r')
hold off
titre='Erreurs relatives des valeurs approchées de x';
title(titre);
disp('done')