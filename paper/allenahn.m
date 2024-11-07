% allencahn.m - solution of Allen-Cahn equation by ETDRK4 scheme
%
%   u_t = 0.01*u_xx + u - u^3 on [-1,1], u(-1)=-1, u(1)=1
%   computation is based on Chebyshev points, so linear term is nondiagonal
%   compare p34.m in Trefethen, "Spectral Methods in MATLAB", SIAM 2000
%   AK Kassam and LN Trefethen, July 2002

% Spatial grid and initial condition:
N = 20;
[D,x] = cheb(N); x = x(2:N);          % spectral differentiation matrix
w = .53*x + .47*sin(-1.5*pi*x) - x; % use w = u-x to make BCs homogeneous
u = [1;w+x;-1];

% Precompute various ETDRK4 matrix quantities:
h = 1/4;                             % time step
M = 32;                              % no. of points for resolvent integral
r = 15*exp(1i*pi*((1:M)-.5)/M);      % points along complex circle
L = D^2; L = .01*L(2:N,2:N);         % 2nd-order differentiation
A = h*L;
E = expm(A); E2 = expm(A/2);
I = eye(N-1); Z = zeros(N-1);
f1 = Z; f2 = Z; f3 = Z; Q = Z;
for j = 1:M
    z = r(j);
    zIA = inv(z*I-A);
    Q = Q + h*zIA*(exp(z/2)-1);
    f1 = f1 + h*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
    f2 = f2 + h*zIA*(2+z+exp(z)*(z-2))/z^2;
    f3 = f3 + h*zIA*(-4-3*z-z^2+exp(z)*(4-z))/z^2;
end
f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);

% Main time-stepping loop:
uu = u; tt = 0;
tmax = 70; nmax = round(tmax/h); nplt = floor((tmax/70)/h);
for n = 1:nmax
    t = n*h;
    Nu = (w+x) - (w+x).^3;
    a = E2*w + Q*Nu;
    Na = a + x - (a+x).^3;
    b = E2*w + Q*Na;
    Nb = b + x - (b+x).^3;
    c = E2*a + Q*(2*Nb-Nu);
    Nc = c + x - (c+x).^3;
    w = E*w + f1*Nu + 2*f2*(Na+Nb) + f3*Nc;
    if mod(n,nplt)==0
	u = [1;w+x;-1];
	uu = [uu,u]; tt = [tt,t];
    end
end

% Plot results:
surf([1;x;-1],tt,uu'), lighting phong, axis tight
view([-45 60]), colormap(cool), light('col',[1 1 0],'pos',[-10 0 10])
