clear all
close all
clc

rng("default")

xstart = -5;
xend = 5;
nplot = 100;

x = linspace(xstart,xend,nplot);

func = @func3;

figure(1)
plot(x, func(x))
xlabel('x-Axis')
ylabel('y-Axis')
title('Function plot')

nstart = 10;
nend   = 400;

count = 0;
for n=nstart:2:nend
    count = count + 1;
    rsmint(count)   = RiemannSumMid(func, xstart, xend, n);
    srint(count)    = Simpson(func, xstart, xend, n);
    mcint(count)    = MonteCarlo(func, xstart, xend, n);
    mcinte(count)   = mean(mcint);

    rsmintEB(count) = RiemannSumMidErr(func, xstart, xend, n);
    srintEB(count)    = SimpsonEB(func, xstart, xend, n);

end

q = integral(func, xstart, xend);
mlbm = ones(size(rsmint)) * q;


figure(2)
hold on
plot(mcint)
plot(mcinte, 'k-d')
plot(rsmint, 'r-x')
plot(srint, 'b-s')
plot(mlbm, 'g-+')

xlabel('Discretization pionts')
ylabel('Integral value')
title('Convergence plot')
legend('MCE','MCEE','RSM','SIM','BM')

figure(3)
semilogy(abs(rsmint - mlbm), 'r-x')
hold on;
semilogy(rsmintEB, 'r:o')

semilogy(abs(srint - mlbm), 'b-s')
semilogy(srintEB,'b:o')

semilogy(abs(mcinte - mlbm), 'k-d')

xlabel('Discretization pionts')
ylabel('Error')
title('Convergence error plot')
legend('RSM-BM','RSM Bound','SIN-BM','SIN Bound', 'MCEE-BM')

function F = MonteCarlo(f, a, b, n)

    x = a+(b-a)*rand(1,n);
    F = (b-a)*sum(f(x))/n;

end 


function bound = SimpsonEB(f,a,b,n)

    % discretization for ddf
    nddf = 5*n*round(b-a) + 500;
    x = linspace(a,b,nddf);
    h = x(2) - x(1);

    % compute second derivative
    q = f(x);
    ddf = (q(1:end-2) - 2*q(2:end-1) + q(3:end))/h^2; 

    % compute M
    M = max(ddf);

    bound = M * (b-a)^5/(180 * n^4);

end

function F = Simpson(f, a, b, n)

    assert(mod(n,2) == 0,'n is not even!')

    % partition
    p = linspace(a,b,n+1);

    % delta x
    dx = p(2) - p(1);

    % eval function
    q = f(p);
    coeff = ones(size(q));
    coeff(2:2:end-1) = 4;
    coeff(3:2:end-2) = 2;

    % Simpson's rule
    F = dx * sum(coeff.*q)/3;

end


function bound = RiemannSumMidErr(f, a, b ,n)

    % discretization for ddf
    nddf = 5*n*round(b-a) + 500;
    x = linspace(a,b,nddf);
    h = x(2) - x(1);

    % compute second derivative
    q = f(x);
    ddf = (q(1:end-2) - 2*q(2:end-1) + q(3:end))/h^2; 

    % compute M
    M = max(ddf);

    bound = M * (b-a)^3/(24 * n^2);

end


function F = RiemannSumMid(f, a, b ,n)

    % partition
    p = linspace(a,b,n+1);

    % Mid points
    x = zeros(n);
    x = (p(1:end-1) + p(2:end))/2;

    % delta x
    dx = p(2) - p(1);

    % function eval
    F = sum(f(x) * dx);

end



function y=func1(x)
    y = x.^2;
end

function y = func2(x)
    y = x.^3 - 3*x + 1;
end

function y = func3(x)
    y = cos(2*x);
end

