function pdf=probability_density_function_of_systematic_error_in_number(n_obs,k,beta,b,systematic_error_probability_density_function)
if nargin<1
    close all
    k=100;
    beta=0.08;
    b=3.2;
    systematic_error_probability_density_function='Gaussian distribution';
end
if strcmp(systematic_error_probability_density_function,'Uniform distribution')
    n_obs_min=(1-beta)^b*k;
    n_obs_max=(1+beta)^b*k;
else
    
    f_x=@(n_obs,k)1/(beta*sqrt(2*pi))*exp(-0.5*(((n_obs./k).^(1/b)-1)/beta).^2);
    n_obs_min=fminbnd(@(n_obs)abs(f_x(n_obs,k)-0.0001),0,k);
    n_obs_max=fminbnd(@(n_obs)abs(f_x(n_obs,k)-0.01),k,k*5,optimset('TolX',1e-1,'Display','off','MaxIter',1e3,'TolFun',1e-12));
%     n_obs_max=fzero(@(n_obs)abs(f_x(n_obs,k)-0.01),k*1.02,optimset('TolX',1e-12,'Display','off','MaxIter',1e6,'TolFun',1e-12));
%     n_obs=round(n_obs_min:n_obs_max);   
end
if nargin<1
    n_obs=n_obs_min:n_obs_max;   
end

if strcmp(systematic_error_probability_density_function,'Uniform distribution')
    f_x=@(k,mu)mu^(1/b)*k.^(1/b-1)/(2*sqrt(3)*beta*b);
    pdf=f_x(n_obs,k);
    index=n_obs<n_obs_min|n_obs>n_obs_max;
    pdf(index)=0;
    
    if k==0
        index=n_obs==0;
        pdf(index)=1;
    end
elseif strcmp(systematic_error_probability_density_function,'Gaussian distribution')
%     f_x=@(n_obs,k)1/(beta*sqrt(2*pi))*exp(-0.5*(((n_obs./k).^(1/b)-1)/beta).^2);
    
%     figure;loglog(n_obs,abs(f_x(n_obs,k)-0.2))
    pdf=f_x(n_obs,k);
    index=n_obs<0;
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
    plot(n_obs,pdf,'r--');
end