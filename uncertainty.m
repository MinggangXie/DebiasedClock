function edf_dsfd=uncertainty(diam,Area,edf_diam)
if nargin<1
    diam=xlsread('data\craters.csv')*1000;%m
    Area=length(diam)/1e-3;
    [diam,~,~,Area]=readdiam('I:\Documents\Paper\2019_Chronology_function\Code2\GUI\data\Giordano_Bruno_All_L.scc');
    diam=diam*1000;
    
    close all;
end
poisson_time=1;
useCSFD=1;
sigma_var=0.05;
if nargin<3
    L=1000;
    DL=logspace(log10(min(diam)*(1-sigma_var*4)),log10(max(diam)*(1+sigma_var*4)),L+1);
    DR=DL(2:end);
    DL=DL(1:end-1);
    edf_diam=sqrt(DL.*DR);
    bin_width=DR-DL;
else
    edf_diam=edf_diam;
    L=length(edf_diam);
    step=edf_diam(2)/edf_diam(1);
    bin_width=edf_diam*(step^0.5-step^-0.5);
end

edf_csfd=zeros(1,L);
edf_dsfd=zeros(1,L);
for i=1:length(diam)
    index=diam(i)-4*sigma_var*diam(i)<edf_diam&edf_diam<diam(i)+4*sigma_var*diam(i);
    y=1/(sqrt(2*pi)*sigma_var*diam(i))*exp(-0.5*((edf_diam(index)-diam(i))/(sigma_var*diam(i))).^2);%/m
    y=y/sum(y.*bin_width(index));
    edf_dsfd(index)=edf_dsfd(index)+y;
    p(i,index)=y;
%     loglog(edf_diam(index),y,'r--');hold on
end
edf_isfd=edf_dsfd.*bin_width;
edf_csfd(end:-1:1)=cumsum(edf_isfd(end:-1:1));
edf_csfd=edf_csfd/edf_csfd(1)*length(diam);
if ~poisson_time
    if useCSFD
        loglog(edf_diam,edf_csfd/Area,'b');hold on
    else
        loglog(edf_diam,edf_dsfd/Area,'b');hold on
    end

    confidence      = 1;
    err=edf_dsfd.^0.5*confidence*edf_csfd(1)/edf_dsfd(1);
    % err=edf_isfd.^0.5./edf_isfd.*edf_csfd*confidence;
    err=edf_isfd.^0.5./bin_width*confidence;
    if useCSFD
        y=[edf_csfd-err, fliplr(edf_csfd+err)]/Area;
    else
        y=[edf_dsfd-err, fliplr(edf_dsfd+err)]/Area;
    end
    yscale=get(gca,'Yscale');
    if strcmp(yscale,'log')
        index=y>0;
    else
        index=y>-inf;
    end
    x=[edf_diam, fliplr(edf_diam)];
    patch(x(index),y(index), 'b', ...
        'EdgeColor', 'None', 'FaceAlpha', 0.1);hold on
    plot(x(index),y(index), 'b--')
end
if nargin<1&&useCSFD
    D=sort(diam);
    N=length(diam):-1:1;
    errorbar(D,N/Area,sqrt(N)/Area,'rs');hold on
end

for i=1:length(edf_diam)
    for j=1:length(diam)
        
    end
end

Dmin=7;%m
N_obs_in=sum(diam>=Dmin);
if Dmin<1000
    N1_km_csfd=logspace(log10(N_obs_in*(Dmin/1000)^3.4/1000),log10(N_obs_in*(Dmin/1000)^3.4*1000),10001);
end

N1_km_width=N1_km_csfd(2:end)-N1_km_csfd(1:end-1);
N1_km_csfd=(N1_km_csfd(2:end).*N1_km_csfd(1:end-1)).^0.5;



[N1_fit,P_in,N1_err_neg,N1_err_pos]=Maximum_likelihood_Michael2016(edf_diam(edf_diam>Dmin)/1000,bin_width(edf_diam>Dmin)/1000,Dmin/1000,N_obs_in,N1_km_csfd,N1_km_width,Area);
% P=P.*probability_number_from_outside(D(D<Dc),0,Dc,sigma_var,N1_km_csfd,Area);
% P_out=probability_number_from_outside(D(D<Dc),N_obs_in-i,Dc,sigma_var,N1_km_csfd,Area);
% P(i+1,:)=P_in.*P_out;
%     N1_fit(i+1,:)=
% -(exp(-(D2 - D)^2/(2*D^2*alpha^2))*(1/D^Dmax - 1/D^Dmin))/log(D)
% figure;
% imagesc(N1_km_csfd,0:N_obs_in,log10(P));colorbar

% plot(edf_diam,(edf_diam/1000).^-3.2*N1_fit,'g-')
[CSFD,DSFD]=getPF(edf_diam/1000,N1_fit);%/km for dsfd


if poisson_time
    N1_fit=density2age(N1_fit,'Neukum et al. (2001)')*1e6;
    N1_km_csfd=density2age(N1_km_csfd,'Neukum et al. (2001)')*1e6;
    N1_err_neg=density2age(N1_err_neg,'Neukum et al. (2001)')*1e6;
    N1_err_pos=density2age(N1_err_pos,'Neukum et al. (2001)')*1e6;
    
    set(gca,'Yscale','log','Xscale','log');

    if useCSFD
        plot(edf_diam,CSFD, 'b--')
        plot(edf_diam(edf_diam>Dmin),CSFD(edf_diam>Dmin), 'b-')
    else
        plot(edf_diam,DSFD/1000, 'b--')
        plot(edf_diam(edf_diam>Dmin),DSFD(edf_diam>Dmin), 'b-')
    end
end
xlabel('Diameter D (m)')
ylabel('Numbers of crater > D (m)')
h=legend(sprintf('$%.1f_{-%.1f}^{+%.1f} ka(1\\sigma)$',N1_fit,N1_fit-N1_err_neg,N1_err_pos-N1_fit));
set(h,'interpreter','latex')

if 1
    h_inset=axes('Position',[0.6,0.3,0.3,0.2]);
%     figure;
    P=P_in;
    if P(end)==1
        index=P>0.001&P<max(P)*0.9999;
    else
        index=P>0.00001;
    end
    semilogx(N1_km_csfd(index),P(index)*100,'r')
    xlim([min(N1_km_csfd(index)) max(N1_km_csfd(index))])
    set(h_inset,'Xscale','log','box','off','ytick',[],'ycolor','w');
end
return
%-------------------------------------------------------------------------------------------
if 0
[N1_fit,P_in,N1_err_neg,N1_err_pos]=Maximum_likelihood_Xie2021(edf_diam(edf_diam>Dmin)/1000,bin_width(edf_diam>Dmin)/1000,Dmin/1000,D(D>Dmin)/1000,N1_km_csfd,N1_km_width,Area);

[CSFD,DSFD]=getPF(edf_diam/1000,N1_fit);%/km for dsfd

if useCSFD
    plot(edf_diam,CSFD, 'b--')
    plot(edf_diam(edf_diam>Dmin),CSFD(edf_diam>Dmin), 'b-')
else
    plot(edf_diam,DSFD/1000, 'b--')
    plot(edf_diam(edf_diam>Dmin),DSFD(edf_diam>Dmin), 'b-')
end
h=legend(sprintf('${%.5f_{-%.5f}}^{+%.5f} ka(1\\sigma)$',N1_fit*1e6,(N1_fit-N1_err_neg)*1e6,(N1_err_pos-N1_fit)*1e6));
set(h,'interpreter','latex')
if 1
    h_inset=axes('Position',[0.6,0.3,0.3,0.2],'box','off');
%     figure;
    P=P_in;
    if P(end)==1
        index=P>0.001&P<max(P)*0.9999;
    else
        index=P>0.00001;
    end
    semilogx(N1_km_csfd(index),P(index)*100,'b--')
    xlim([min(N1_km_csfd(index)) max(N1_km_csfd(index))])
    set(h_inset,'Xscale','log');
end
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function [N1_fit,P,N1_err_neg,N1_err_pos]=Maximum_likelihood_Xie2021_(edf_diam,bin_width,D_min,D,N1_km_csfd,N1_km_width,Area)%???
% C=exp(-Area*N1_km_csfd.*(D_min/1000)^-3.2);
% for i=1:length(C)
%     if C(i)<1e-10
%         P(i)=Area*N1_km_csfd(i).*(D_min/1000)^-3.2*N1_km_csfd(i).^Ncum_edf;
%     else
%         P(i)=C(i)*N1_km_csfd(i).^Ncum_edf;
%     end
% end
N_obs_in=length(D);
if 1
    N_Dmin=Area*getPF(D_min,N1_km_csfd);%number, not density
else
    CSFD=@(D,N1_km_csfd)N1_km_csfd*D.^-3.2;%km
    lambda=Area*CSFD(D_min,N1_km_csfd);
end
edf_diam_L=edf_diam*(edf_diam(1)/edf_diam(2))^0.5;
index_D=zeros(size(edf_diam));
for i=1:length(D)
    index=find(D(i)>=edf_diam_L,1,'last');
    index_D(index)=1;
end
index_D_comp=~index_D;%complement set
index_D=~index_D_comp;

alpha=0.1;
for i=1:length(N1_km_csfd)
    P_poisson(i)=poisspdf(N_obs_in,N_Dmin(i));
    
    [CSFD,DSFD]=getPF(edf_diam,N1_km_csfd(i));
    c=sum(log(DSFD./CSFD))+sum(log(bin_width));
    if P_poisson(i)==0
        P_log(i)=-inf;
    else
        P_log(i)=c+log(P_poisson(i));
    end
    lambda=Area*DSFD.*bin_width;
    err_gaussian=sum(log(1./(sqrt(2*pi)*DSFD.*edf_diam)))-sum(lambda(index_D_comp).^2./(2*(alpha*edf_diam(index_D_comp).*DSFD(index_D_comp))).^2)...
    -sum((1-lambda(index_D)).^2./(2*(alpha*edf_diam(index_D).*DSFD(index_D))).^2);
    P_log(i)=P_log(i)+err_gaussian;
%     P_log(i)=N_obs_in*log(Area)+sum(log(bin_width))-N_obs_in+sum(log(DSFD));
end
P_log=P_log-min(P_log(P_log>-inf));
P=exp(P_log);

P_cum=cumsum(P.*N1_km_width);
P=P/P_cum(end);
P_cum=P_cum/P_cum(end);
[P_cum2,index]=unique(P_cum);
N1_fit=interp1(P_cum2,N1_km_csfd(index),0.5);
N1_err_neg=interp1(P_cum2,N1_km_csfd(index),(1-0.68)/2);
N1_err_pos=interp1(P_cum2,N1_km_csfd(index),1-(1-0.68)/2);


function uncertianty_new(edf_diam,edf_N,N1_km_csfd)
PROB=edf_N/sum(edf_N);
p=mnpdf(edf_N,PROB);
plot(edf_diam,p);

function number_uncertainty(D,N)
if nargin<1
    D=xlsread('data\craters.csv')*1000;
end
% p=


function P=probability_number_from_outside(D,N,Dc,sigma_var,N1_km_csfd,Area)
if nargin<5
    N_expected=sum(erfc((D-Dc)./(sqrt(2)*sigma_var*D)));
    P=poisspdf(N,N_expected);
else
    DL=logspace(log10(Dc/2),log10(Dc),101)/1000;%km
    DR=DL(2:end);
    DL=DL(1:end-1);
    D=sqrt(DL.*DR);
    for i=1:length(N1_km_csfd)
        [~,DSFD]=getPF(D,N1_km_csfd(i));
        N_expected=sum(Area*DSFD.*(DR-DL).*erfc((D-Dc)./(sqrt(2)*sigma_var*D))/2);
        P(i)=poisspdf(N,N_expected);
    end
end


function [CSFD,DSFD,ISFD]=getPF(D_km,N1_km_csfd)%D,km;
CSFD=N1_km_csfd*D_km.^-3.2;%km
DSFD=3.2*N1_km_csfd*D_km.^-4.2;%km













