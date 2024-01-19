close all 
clear all
clc

rng("default")

% Looking at the function to be integrated
n      = 1000;
xstart = -2;
xend   = 2;
x      = linspace(xstart,xend,n);

func = @integrant1;

figure(1)
plot(x,func(x))
xlabel('x-axis')
ylabel('y-axis')
title('Function plot')

% Starting numerical integration 

nstart = round(xend - xstart);
nend   = 1000*nstart;

count = 0;
for n = nstart:2:nend

    count = count + 1;     % Riemann Sum mid point
    RSM(count) = RiemannMid(func, xstart, xend, n);
    RSME(count) = RiemannMidEB(func, xstart, xend, n);

    % Simpson rule
    Sim(count) = Simpson(func, xstart, xend, n);
    SimE(count) = SimpsonEB(func, xstart, xend, n);

    % Monte Carlo
    MCint(count) = MC(func, xstart, xend, n);

end

% MATLAB benchmark
MLint = ones(size(RSM)) * integral(func, xstart, xend);

figure(2)
hold on
plot(MCint, "g-s")
plot(MLint, "k-|" )
plot(RSM, "r-+" )
plot(Sim, "b-x")
legend('MC','MATLAB','RMS', 'Sim')

figure(3)
hold on
plot(abs(MLint - RSM), "r-+")
plot(RSME,"r:o")

plot(abs(MLint - Sim), "b-+")
plot(SimE,"b:o")

plot(abs(MLint - MCint), "g-+")
legend('BM - RSM','RMS', 'BM - Sim', 'Sim', 'BM - MC')

function y = integrant1(x)
    y = x.^2;
end

function y = integrant2(x)
    y = x.^3 - 3*x + 1;
end

function y = integrant3(x)
    y = cos(2*x);
end

function y = integrant4(x)
    y = exp(-x.^2).*cos(5*x);
end

function y = integrant(x)
    y = sin(x).^2;
end

function F = MC(f, a, b, n)

    % Uniformly distributed RV on [a,b]
    x = a + (b-a)*rand(1,n);

    % Monte-Carlo Estimator  
    F = (b - a) * sum(f(x)) / n;

end


function F = Simpson(f, a, b, n)

    assert(mod(n,2) == 0)

    % Generate partition
    p = linspace(a,b,n+1);

    % compute delta x
    dx = (b-a)/n ;

    % Function values at partition points
    y = f(p);

    % Simpson's coefficients
    coefficients = ones(1, n + 1);
    coefficients(2:2:end-1) = 4;
    coefficients(3:2:end-2) = 2;

    % Integral approximation
    F = dx/3 * sum(y .* coefficients);

end

function F = RiemannMid(f, a, b, n)

    % Generate partition
    p = linspace(a,b,n+1);
    
    % Compute mid points
    x = zeros(1, n);
    x = (p(1:end-1) + p(2:end))/2;

    % compute delta x
    dx = (b-a)/n ;

    F = sum(f(x).* dx);

end

function bound = SimpsonEB(f, a, b, n)

    x = linspace(a,b, 2 * n + 1);
    h = x(2)- x(1);

    % finite difference second order
    ddf = zeros(1, length(x)-2);
    ddf = (f(x(3:end)) - 2*f(x(2:end-1)) + f(x(1:end-2)))/h^2;

    M = max(ddf);

    bound = M * (b -a)^5 / (180 * n^4);

end
    
function bound = RiemannMidEB(f, a, b, n)

    x = linspace(a,b, 2 * n + 1);
    h = x(2)- x(1);

    % finite difference second order
    ddf = zeros(1, n-1);
    ddf = (f(x(3:end)) - 2*f(x(2:end-1)) + f(x(1:end-2)))/h^2;

    M = max(ddf);

    bound = M * (b -a)^3 / (24 * n^2);

end
