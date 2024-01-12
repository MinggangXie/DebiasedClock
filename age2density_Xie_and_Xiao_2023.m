function N1_debiased=age2density_Xie_and_Xiao_2023(t)
if nargin<1
    t=linspace(0,4.5,1000);
end
N1_debiased=7.829e-4*t+51.11*(0.04298*exp((t-4.5)/6.208e-3)+0.5766*exp((t-4.5)/2.707e-2)+0.3779*exp((t-4.5)/7.933e-2));

if nargin<1
    figure;
    semilogy(t,N1_debiased);
end