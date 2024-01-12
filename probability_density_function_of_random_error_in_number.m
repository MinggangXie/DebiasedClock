function pdf=probability_density_function_of_random_error_in_number(n_obs,k,sigma_alpha,b,systematic_error_probability_density_function)
if nargin<1
    close all
    k=10;
    sigma_alpha=0.02;
    b=3.2;
    systematic_error_probability_density_function='Gaussian distribution';
end
if strcmp(systematic_error_probability_density_function,'Uniform distribution')
    beep;
    fprintf(2,'The code for Uniform distribution is not ready');
    return;
end
f_nobs_expected=@(b,alpha)(1+0.25*b^2.3*alpha^2)*k;
if strcmp(systematic_error_probability_density_function,'Uniform distribution')
    n_obs_min=(1-sigma_alpha)^b*k;
    n_obs_max=(1+sigma_alpha)^b*k;
else
    n_obs_min=f_nobs_expected(b,sigma_alpha)/5;
    n_obs_max=f_nobs_expected(b,sigma_alpha)*5;
end
if nargin<1
    n_obs=round(0:n_obs_max);   
end

if strcmp(systematic_error_probability_density_function,'Uniform distribution')
    f_x=@(k,mu)mu^(1/b)*k.^(1/b-1)/(2*sqrt(3)*sigma_alpha*b);
    pdf=f_x(n_obs,k);
    index=n_obs<n_obs_min|n_obs>n_obs_max;
    pdf(index)=0;
    
    if k==0
        index=n_obs==0;
        pdf(index)=1;
    end
elseif strcmp(systematic_error_probability_density_function,'Gaussian distribution')
    f_x=@(n_obs,k)1/(sigma_alpha*sqrt(2*pi))*exp(-0.5*(((n_obs./k-1).^0.5/(0.5*b^1.15))/sigma_alpha).^2);
    
%     n_obs_min=k;
%     n_obs_max=(1+0.25*b^2.3*(sigma_alpha*3.5)^2)*k;
%     n_obs=linspace(n_obs_min,n_obs_max,20);   
    
    pdf=f_x(n_obs,k);
    index=n_obs<k;
    pdf(index)=0;
    
    if k==0
        index=n_obs==0;
        pdf(index)=1;
    end
end
if sum(pdf)>0
    pdf=pdf/sum(pdf);
end
if nargin<1
%     pdf_poiss=poisspdf(n_obs, k);
%     plot(n_obs,pdf_poiss,'k-');hold on
    plot(n_obs,pdf,'r+--');
end