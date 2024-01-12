function [lower_bound_1sigma,upper_bound_1sigma]=get_uncertainty(confidence_interval,x,pdf)
% Parameters
% confidence_interval = 0.95; Confidence interval

% Calculate the uncertainty
ecdf_values = cumsum(pdf) / sum(pdf); % Empirical cumulative distribution function (CDF)
lower_percentile = (1 - confidence_interval) / 2; % Lower percentile for the confidence interval
upper_percentile = 1 - lower_percentile; % Upper percentile for the confidence interval
[ecdf_values_unique,idx]=unique(ecdf_values);
x_unique=x(idx);
lower_bound_1sigma = interp1(ecdf_values_unique, x_unique, lower_percentile, 'linear', 'extrap'); % bound of confidence interval
if lower_bound_1sigma<0
    lower_bound_1sigma=0;
end
upper_bound_1sigma = interp1(ecdf_values_unique, x_unique, upper_percentile, 'linear', 'extrap'); % Upper bound of confidence interval

% uncertainty = upper_bound_1sigma - lower_bound_1sigma; % Uncertainty calculation

% Display the uncertainty
% disp(['The uncertainty of the value within the confidence interval is: ' num2str(ertainty)]);

