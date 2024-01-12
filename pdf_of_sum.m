function [x,pdf_PSSR]=pdf_of_sum(lambda,beta,b,systematic_error_probability_density_function)
if nargin<1
    lambda=10;
    beta=0.08;
    b=3.2;
    systematic_error_probability_density_function='Gaussian distribution';
end
dx=1;
x = floor(lambda/20):dx:lambda*5; % Range of values for the sum

%----------------------------------------------------------Poisson+systematic error ----------------------------------------------------------------------
pdf_PSS = zeros(size(x));
for k = 0:max(x)
    prob_poisson = poisspdf(k, lambda); % Probability of Poisson variable taking value k
    %     pdf_gaussian = normpdf(x-k, mu, sigma); % pdf_PSS of Gaussian variable at x-k
    pdf_gaussian=probability_density_function_of_systematic_error_in_number(x,k,beta,b,systematic_error_probability_density_function);%
    %     plot(x,pdf_PSS)
    pdf_PSS = pdf_PSS + prob_poisson * pdf_gaussian; % Add weighted probabilities
end
% Normalize the pdf_PSS
pdf_PSS = pdf_PSS / sum(pdf_PSS);

if nargin<1
    plot(x,pdf_PSS,'k-');hold on
end
%---------------------------------------------------------------(Poisson+systematic error)+random error-----------------------------------------------------------------
if 1 %random error is neglected
    pdf_PSSR=pdf_PSS;
else
    pdf_PSSR=pdf_PSS;
    for k = 0:max(x)
        pdf_random_error=probability_density_function_of_random_error_in_number(x,k,beta,b,systematic_error_probability_density_function);%
        %     plot(x,pdf_SS)
        pdf_PSSR = pdf_PSSR + pdf_PSS(k+1) * pdf_random_error; % Add weighted probabilities
    end
    % Normalize the pdf_PSS
    pdf_PSSR = pdf_PSSR / sum(pdf_PSSR);
end
if nargin<1
    plot(x,pdf_PSSR,'r--')
end