close all
clear all
clc
rng("default")

dim = 100;
its = 2000;
L = Laplacian_2D(dim, 1);

%
tic
B = Arnoldi(L, 10);
disp('Elapsed time Arnoldi:')
disp(toc)

% Solving EV problem with QR "quick and dirty" (RR)
tic
M = B.'* L *B;
%surf(M);view(0,90)
e_vals = qr_eigs(M,its);
disp('Elapsed time RR:')
disp(toc)

% Solving EV problem with sketched QR "quick and dirty" (sRR)
tic
e_vals_srr = qr_eigs_sketch(L, B, 1, its);
disp('Elapsed time s-RR:')
disp(toc)

% MATLAB eig
target_vals = sort(eig(M));

disp('### eigenvalue differences ###')
disp('Rayleigh-Ritz with QR:')
disp([target_vals, e_vals,abs(target_vals - e_vals)])
disp('Sketched Rayleigh-Ritz with QR:')
disp([target_vals, e_vals_srr, abs(target_vals - e_vals_srr)])

hold on
plot(target_vals)
plot(e_vals,'r--')
plot(e_vals_srr, 'kx')

function [evals] = qr_eigs_sketch(A, B, type, its)
%tol = 1e-6;

m = size(A,1);
d = size(B,2);
s = 4*d;

% sketching B and A*B
if type == 1
    % Gauss
    S = rand(s,m);   
    SB = S*B;
    SAB = S*A*B;
elseif type == 2
    % SRTT
    [SB, SAB] = SRTT(A,B,s);
elseif type == 3
    % SEE
    [SB, SAB] = SSE(A,B,s);
else
    disp('Undefined sketch')
end

[U,T] = qr(SB);
T_inv = inv_T(T);

%M = pinv(SB) * SAB;
%M = pinv(T) * U.' * SAB;
M = T_inv * U.' *SAB;

diff = inf;
val = inf;
%while diff > tol
for i=1:its
    [Q,R] = qr(M);
    M = R * Q;
    evals = sort(diag(M));
    val_new = max(abs(evals));
    diff = abs(val_new - val);
    val = val_new;
end

end

function U_inv = inv_T(U)
    [m, n] = size(U);  
    U_inv = zeros(n, m);
    
    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1; 
        U_inv(:, i) = backward_substitution(U, e);
    end
end

function x = backward_substitution(U, b)
    n = size(U,2);
    x = zeros(n, 1);
 
    for i = n:-1:1
        x(i) = (b(i) - U(i, i+1:end) * x(i+1:end)) / U(i, i);
    end
end

function [SB,SAB] = SSE(A,B,d)

zeta = 8;

m = size(B,1);
S = zeros(d,m);

for i = 1:m
    idx = randsample(d,zeta);
    signs = 2 * randi(2,zeta,1) - 3;
    S(idx,i) = signs/sqrt(zeta);
end

SB = S * B;
SAB = S * A * B;
end

function [SB,SAB] = SRTT(A,B,d)
% Sketching B, and A*B

[n,m] = size(B);

signs = 2 * randi(2,n,1) - 3;
idx = randsample(n,d);

SB = signs.* B;
SB = dct(SB);
SB = SB(idx,:);
SB = sqrt(n/d) * SB;

AB = A*B;
cAB = signs.* AB;
cAB = dct(cAB);
cAB = cAB(idx,:);
cAB = sqrt(n/d) * cAB;

SAB = cAB;

end

function [evals] = qr_eigs(A, its)

%tol = 1e-6;

diff = inf;
val = inf;
%while diff > tol
for i=1:its
    [Q,R] = qr(A);
    A = R * Q;
    evals = sort(diag(A));
    val_new = max(abs(evals));
    diff = abs(val_new - val);
    val = val_new;
end
end

function [L] = Laplacian_2D(N, dx)

    L_1D = -2*eye(N); 
    L_1D = L_1D + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    L_1D = L_1D / dx^2;

    L = kron(eye(N), L_1D) + kron(L_1D, eye(N));
end

function [Q] = Arnoldi(A,p)

    r = randn(size(A,1),1);
    q = r/norm(r);
    Q = q;
    for i = 1:p-1
        Aq = A*q;
        w = Aq - Q*(Q.'*Aq);
        q = w/norm(w);
        Q = [Q,q];
    end    

end