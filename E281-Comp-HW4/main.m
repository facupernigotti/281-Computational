


%% Problem 2
% Define the range of values for r
clear all ; close all ; clc ;

tic ;
r_values = 0:0.001:0.04;

% Initialize an array to store the results
results = zeros(size(r_values));

% Calculate the results for each value of r
for i = 1:numel(r_values)
    results(i) = av_wealth_calc(r_values(i));
end

% Plot the results
plot(r_values, results);
xlabel('r');
ylabel('Average Wealth');
title('Average Wealth vs. r');

toc ;
