function []=test_EmpiricalPDF()
clear all
close all 
clc

function y = true_pdf(x)
    y = 1 / sqrt(2 * pi) * exp(-0.5 * x.^2);
end

function output = empirical_pdf(x,samples, num_bins)
    min_val = min(x);
    max_val = max(x);
    edges = linspace(min_val, max_val, num_bins+1);
    diff = edges(2)- edges(1);

    counts = zeros(1, num_bins);    
    for i = 1:num_bins
        counts(i) = sum(samples >= edges(i) & samples < edges(i+1))/(size(samples,1) * diff);
    end

    num_values = size(x,2);
    output = zeros(1, num_values); 
    for i =1:num_values
        for j =1:num_bins
            if (x(i) >= edges(j) & x(i) < edges(j+1))
                output(i) = counts(j);
            end
        end
    end
end

% Sample the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw random data
rng('default'); % For reproducibility
num_samples = 600;
samples = randn(num_samples, 1); % Standard normal distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse transform sampling
data = randn(num_samples, 1); % Sample data from an unknown distribution

% Sort the data
sorted_data = sort(data);

% Estimate the empirical CDF
empirical_cdf = (1:length(sorted_data)) / length(sorted_data);

% Generate uniform random numbers
u = rand(num_samples, 1);

% Use the empirical CDF to transform the uniform random numbers into samples
samples2 = interp1(empirical_cdf, sorted_data, u, 'linear', 'extrap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_vals = 1000;
x_values = linspace(-3, 3, num_vals);
true_pdf_values = true_pdf(x_values);
empirical_pdf_values = empirical_pdf(x_values, samples, 20);
empirical_pdf_values2 = empirical_pdf(x_values, samples2, 20);

% Plot the true PDF and empirical PDF
plot(x_values, true_pdf_values, 'b-')
hold on;
plot(x_values, empirical_pdf_values, 'r--');
plot(x_values, empirical_pdf_values2, 'k--');


xlabel('x');
ylabel('Probability Density');
title('Comparison of True and Empirical Probability Density Functions');
legend('True PDF', 'Rand Sample', 'IT Sample', 'Location', 'best');
grid on;

end

