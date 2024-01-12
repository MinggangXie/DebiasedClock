function t=density2age(N1,ChronologyFunction)
%N1, km^{-2}
%t,Ga
if nargin<1
    N1=0.003;
    ChronologyFunction='Xie and Xiao (2023)';
end
if N1==0
    t=0;
    return;
end
t=zeros(size(N1));
switch ChronologyFunction
    case 'Xie and Xiao (2023)'
        t=fminbnd(@(t)abs(age2density_Xie_and_Xiao_2023(t)-N1),0,4.5,optimset('TolX',N1/7.829e-4*1e-5,'Display','off'));
    case {'Neukum (1983)','Neukum et al. (2001)'}
        t=(4611686018427387904*N1)/3864592883442151 - (100*lambertw(0, (21334559601841101*exp((114139228956077850624*N1)/13802117440864825 ...
        + 21334559601841101/47423714419252509865410560))/47423714419252509865410560))/693 + 153929001456285/2371185720962625493270528;
    case 'Xie et al. (2018) (Default)'
        for i=1:length(N1)
            f=@(x)abs(1.0033e-11*(exp(5.5757*x)-1)+9.986e-05*x+0.00010035*x.^2+(0.00017336*(erf(0.12052/0.061949)-erf((0.12052-x)/0.061949)))/2+(0.00036988*(erf(0.45/0.06)-erf((0.45-x)/0.06)))/2-N1(i));
        %     f=@(x)abs(1.3077e-11*(exp(5.4761*x)-1)+0.0003016*x+0.00016622*x.^2-N1(i));
            t(i)=minHJ(f,0,5,1e-10);
        end
    otherwise 
        b=str2num(PF_type(3:end));
        str_productionfunction=['a0*x^-' PF_type(3:end)];
end

