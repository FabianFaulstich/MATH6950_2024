close all
clear all
clc
rng("default")

dim = 100;
its = 2000;
p = 20;
L = Laplace(dim,1);
L = sparse(L);

%spec_L = sort(eig(L));
%plot(spec_L)

% Computing Arnoldi yields upper Hessenberg
tic
B = Arnoldi(L, p);
disp('Elapsed time for Arnoldi:')
disp(toc)

% Compute the specrtrum using QR "quick and dirty"
tic
M = B.' * L  * B;
evals_qr = eigs_qr(M, its);
disp('Elapsed time for QR:')
disp(toc)

% Compute the spectrum using s-RR
tic
evals_srr = eigs_srr(L, B, its);
disp('Elapsed time for s-RR')
disp(toc)

% Compare spectrum of M
evals = sort(eig(M));
disp('Coparison of spectra using QR:')
disp([evals, evals_qr, abs(evals - evals_qr)])
disp('Coparison of spectra using s-RR:')
disp([evals, evals_srr, abs(evals - evals_srr)])

hold on
plot(evals)
plot(evals_qr, 'r--')
plot(evals_srr, 'b*')

function [evals] = eigs_srr(A, B, its)
    
    [m,p]= size(B);
    d = 4*p;

    % sketching
    %S = rand(d,m);
    %SB = S*B;
    %SAB = S*A*B;
    [SB, SAB] = SRTT(A, B, d);
    %[SB, SAB] = SSE(A, B, d);

    M = pinv(SB)*SAB;

    for i=1:its
        [Q,R] = qr(M);
        M = R * Q;
    end
    evals = sort(diag(M));    
end




function [SB, SAB] = SSE(A,B,p)

    m = size(A,1);
    zeta = 8;
    
    S = spalloc(p,m,zeta*m);
    
    for i=1:m
        idx = randsample(p,zeta);
        signs = 2*randi(2,zeta,1)-3;
        for j=1:zeta
            S(idx(j),i) = signs(j)/sqrt(zeta);
        end
    end

    SB = S*B;
    SAB = S*A;
    SAB = SAB * B;

end

function [SB, SAB] = SRTT(A,B,d)
    
    n = size(B,1);
    AB = A*B;

    signs = 2*randi(2,n,1)-3;
    idx = randsample(n,d);

    SB = signs.*B;
    SB = dct(SB);
    SB = SB(idx,:)/sqrt(n/d);

    SAB = signs.*AB;
    SAB = dct(SAB);
    SAB = SAB(idx,:)/sqrt(n/d);

end

function [evals] = eigs_qr(A, its)
    
for i=1:its
   [Q,R] = qr(A);
   A = R * Q;
end

evals = sort(diag(A));
end

function [B] = Arnoldi(A, p)

n = size(A,2);

r = rand(n,1);
b = r/norm(r);

B = b; 

for i = 1:p-1
    w = A*b;
    b = w - B*(B.'*w);
    b = b/norm(b);
    B = [B,b];
end

end


function [L] = Laplace(N,dx)

    L_1D = -2*eye(N); 
    L_1D = L_1D + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    L_1D = L_1D / dx^2;

    L = kron(eye(N), L_1D) + kron(L_1D, eye(N));

end