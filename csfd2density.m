function [N1_fitting,N1_err,DL_fit,cumN_fit]=csfd2density(DL,cumN,err,area,PF_type)
%[age,DL_fit,cumN_fit]=csfd2age(DL,cumN,err,area)
if nargin<1
    close all
    folder='I:\Documents\Paper\Small craters\ArcGIS\Serenitatis\划分地层V1\shp\Area_i\';
    filename_area=[folder  '4.dbf'];
    data=shaperead(filename_area);
    area=data.Area;
    filename=[folder 'Craters_4.dbf'];
    [DL,cumN,err]=shp2cum(filename,area,'pseudo-log');%2^(1/8)
    
    PF_type='Neukum et al. (2001)';
    PF_type='Xie et al. (2019) (Default)';
%     PF_type='b=3';%N=D^-b
end
for i=length(DL):-1:2
    if cumN(i)==cumN(i-1)
        cumN(i-1)=[];
        DL(i-1)=[];
        err(i-1)=[];
    end
end
if nargin<1
    errorbarloglog(DL,cumN,err,'go');hold on
end
switch PF_type
    case 'Neukum et al. (2001)'
        if DL(1)<0.01
            fprintf(2,'Warning: DL<0.01 km\n');
            index=DL>0.01;
            DL=DL(index);
            cumN=cumN(index);
            err=err(index);
        end
end
% errorbarloglog(DL,cumN,err,'ro');hold on
if length(DL(:,1))==1
    DL=DL';
end
if length(cumN(:,1))==1
    cumN=cumN';
end
if length(err(:,1))==1
    err=err';
end
switch PF_type
    case 'Neukum (1983)'
        str_productionfunction='a0-3.6269*x.^1+0.43662*x.^2+0.79347*x.^3+0.086468*x.^4-0.26485*x.^5-0.066382*x.^6+0.037923*x.^7+0.010596*x.^8-0.0022496*x.^9-0.00051797*x.^10+0.0000397*x.^11';
    case 'Neukum et al. (2001)'
        str_productionfunction='a0-3.557528*x.^1+0.781027*x.^2+1.021521*x.^3-0.156012*x.^4-0.444058*x.^5+0.019977*x.^6+0.086850*x.^7-0.005874*x.^8-0.006809*x.^9+0.000825*x.^10+5.54e-05*x.^11';
    case 'Xie_2017'
        str_productionfunction='a0-3.2*x';
    case 'Xie et al. (2019) (Default)'
        str_productionfunction='a0-3.2*x';%y=ax^-b,logy=a-blog10(x)
%         str_productionfunction='a0-0.022962*x.^9-0.21762*x.^8-0.90205*x.^7-2.1586*x.^6-3.121*x.^5-2.0581*x.^4+0.98224*x.^3+2.3706*x.^2-2.5938*x';
%         data=load('data\dataset.mat');
    otherwise 
        b=str2num(PF_type(3:end));
        str_productionfunction=['a0*x^-' PF_type(3:end)];
end
switch PF_type
    case 'Xie et al. (2019) (Default)'
        [N1_fitting,N1_err,DL_fit,cumN_fit]=Chi_square_fitting(DL',cumN',err');
        
        if 0
        opts=fitoptions('Method','NonlinearLeastSquares');%设置拟合用的方法和初始值
        opts.StartPoint=log10(cumN(1));
        opts.TolFun=1e-10;%default
        opts.Weights = (err./(cumN*log(10))).^-2;
        ft = fittype(str_productionfunction,'options',opts);%设置拟合函数为高斯函数

        cfun = fit(log10(DL),log10(cumN),ft);%拟合
        N01=10^cfun.a0;
        index=floor(log10(N01/data.N01)/log10(data.Dstep))+1;
        OPF=interp1(data.DL,data.OPF(index,:),DL);
        N1=OPF./cumN
        
        DL_fit=logspace(log10(DL(1)),3,1e4);
        cumN_fit=f(log10(DL_fit),a0);
        end
    otherwise     
        a0_guess=log10(interp1(DL,cumN,1,'pchip','extrap'));
        index=find(DL<1,1,'last');
        index=1;
        a0_guess=log10(cumN(index)*(DL(index))^3.2);
        % a0=a0_guess;
        opts=fitoptions('Method','NonlinearLeastSquares');%设置拟合用的方法和初始值
        opts.StartPoint=a0_guess;
        opts.MaxFunEvals=4e3;
        opts.MaxIter=2e3;
        opts.TolFun=1e-25;
        opts.TolX=1e-25;
        if 1
            opts.Weights = err.^-2;
            ft = fittype(sprintf('10^(%s)',str_productionfunction),'options',opts);%设置拟合函数为高斯函数

            [cfun,gof] = fit(log10(DL),cumN,ft);%拟合
            
            a0=cfun.a0;

            if abs(a0_guess-a0)<1e-6
                fprintf(2,'a0 may be the same as a0_guess\n');
            end
            c=confint(cfun,0.6827);
            N1_fitting=10^a0;
            f=str2func(['@(x,a0)10.^(',str_productionfunction ')']);
            N1_err=max(N1_fitting-10^c(1),10^c(2)-N1_fitting);
        else
            opts.Weights = (err./(cumN*log(10))).^-2;
            ft = fittype(str_productionfunction,'options',opts);%设置拟合函数为高斯函数

            [cfun,gof] = fit(log10(DL),log10(cumN),ft);%拟合
            
            a0=cfun.a0;

            if abs(a0_guess-a0)<1e-6
                fprintf(2,'a0 may be the same as a0_guess\n');
            end

            N1_fitting=10^a0;
            f=str2func(['@(x,a0)10.^(',str_productionfunction ')']);
            N1_err=err(1)/f(log10(DL(1)),a0)*f(0,a0);
        end
        
end
switch PF_type
    case {'Neukum (1983)','Neukum et al. (2001)'}
        DL_fit=logspace(log10(0.01),log10(300),1e4);
%             cumN_fit=productionfunction(a0,DL_fit);
        cumN_fit=f(log10(DL_fit),a0);
    case 'Xie_2017'
        DL_fit=logspace(log10(DL(1)),3,1e4);
        cumN_fit=10.^(log10(a0)-3.2*log10(DL_fit));
    case 'Xie et al. (2019) (Default)'
%         DL_fit=logspace(log10(DL(1)),3,1e4);
%         cumN_fit=f(log10(DL_fit),a0);
    otherwise 
        DL_fit=logspace(log10(DL(1)),3,1e3);
        cumN_fit=a0*DL_fit.^-b;
end

% N1_measured_err=(N1_fitting*area)^0.5/area;
% age_lower=density2age(N1_fitting-N1_measured_err);
% age_upper=density2age(N1_fitting+N1_measured_err);
% 
% age(2)=age_lower;
% age(3)=age_upper;
if nargout>1||nargin<1
%         cumN_fit=double(subs(productionfunction,DL_fit));
    if nargin<1
        loglog(DL_fit,cumN_fit,'r-');
    end
end
% productionfunction=@(a,x)10.^(a-3.5575*log10(x).^1+0.78103*log10(x).^2+1.0215*log10(x).^3-0.15601*log10(x).^4-...
%         0.44406*log10(x).^5+0.019977*log10(x).^6+0.08685*log10(x).^7-0.005874*log10(x).^8-...
%         0.006809*log10(x).^9+0.000825*log10(x).^10+5.54e-05*log10(x).^11);
% 
%     a0=log10(interp1(DL,cumN,1));
%     a0= lsqcurvefit(productionfunction,a0,DL,cumN,-100,100);%,'Weight', (area*err').^-2
%     age=density2age(10^a0);
%     
end
% 
% 
function f=productionfunction(a0,x)
    f=10.^(a0-3.5575*log10(x).^1+0.78103*log10(x).^2+1.0215*log10(x).^3-0.15601*log10(x).^4-...
        0.44406*log10(x).^5+0.019977*log10(x).^6+0.08685*log10(x).^7-0.005874*log10(x).^8-...
        0.006809*log10(x).^9+0.000825*log10(x).^10+5.54e-05*log10(x).^11);
%     syms x
%     aj=[-3.0876 -3.557528 0.781027 1.021521 -0.156012 -0.444058 0.019977 0.086850 -0.005874 -0.006809 8.25e-4 5.54e-5];
%     str1=[];
%     for i=2:length(aj)
%         if aj(i)>0
%             str1=[str1,'+',num2str(aj(i)),'*log10(x).^', num2str(i-1)];%NPF=NPF+aj(i).*(log10(DL)).^(i-1);
%         else
%             str1=[str1,num2str(aj(i)),'*log10(x).^', num2str(i-1)];%NPF=NPF+aj(i).*(log10(DL)).^(i-1);
%         end
%     end
%     f=['10.^(a0',str1,')'];
end

function [N1_fitting,N1_err,DL_fit,cumN_fit]=Chi_square_fitting(DL,N,err_N,N_PF,N1_PF,DL_PF,PF_slope)
if nargin<4
    N1_PF=logspace(-10,-1,1000);
    DL_PF=logspace(-3,3,500);%km
    
    for i=1:length(N1_PF)
        N_PF(i,:)=N1_PF(i)*DL_PF.^-3.4;
        PF_slope(i,:)=-3.4*N1_PF(i)*DL_PF.^-4.4;
    end
    
end
debug=0;

bin_width=diff(DL);
bin_width=[bin_width,bin_width(end)];%??????????????????

N_Dmin=interp1(DL_PF,N_PF(1,:),DL(1));
PF=interp1(DL_PF,N_PF(1,:),DL);
N_est=sum(PF.*N)/sum(PF.^2)*N1_PF(1);%derivation given in N_est.wmf in Manuscript folder.
err_N1_fit_est=(N1_PF(1)/N_Dmin)*err_N(1);


step=N1_PF(2)/N1_PF(1);
if N_est<N1_PF(1)
    fprintf(2,'N_est<min(N1_PF)\n\n');
    beep
    return;
elseif N_est>N1_PF(end)
    fprintf(2,'N_est>max(N1_PF)\n\n');
    beep
    return;
else
    index_est=log(N_est/N1_PF(1))/log(step)+1;%N_est=N1_PF(1)*step^(k-1) => k=log(N_est/N1_PF(1))/log(step)+1
    
    N_estL=N_est-50*err_N1_fit_est;
    if N_estL<N1_PF(1)
        index_L=1;
    else
        index_L=log(N_estL/N1_PF(1))/log(step)+1;
    end
    N_estU=N_est+1000*err_N1_fit_est;
    if N_estU>N1_PF(end)
        index_U=length(N1_PF);
    else
        index_U=log(N_estU/N1_PF(1))/log(step)+1;
    end
    index_N1=round(index_L:index_U);
end
if isempty(index_N1)
    index_N1=round(max(1,index_est-5):min(index_est+5,length(N1_PF)));
end
if debug
    figure;
end
sigma=0.1*DL;
sigma=0;
for i=1:length(index_N1)
    PF=interp1(DL_PF,N_PF(index_N1(i),:),DL);
    slope=interp1(DL_PF,PF_slope(index_N1(i),:),DL,'next');
%     weight=err_N.^2+slope.^2.*sigma.^2;%Barker 1974
    weight=err_N.^2*1.5;%Barker 1974
    chi_square(i)=sum((PF-N).^2/(weight));
    
    if debug
        loglog(DL_PF,N_PF(index_N1(i),:));hold on
    end
end


[min_chi_square,index_chi_square]=min(chi_square);
N1_fitting=N1_PF(index_N1(index_chi_square));

index=find(chi_square<min_chi_square+2.3,1,'first');
err_L=N1_PF(index_N1(index));
index=find(chi_square<min_chi_square+2.3,1,'last');
err_U=N1_PF(index_N1(index));
N1_err=(err_U-err_L)/2;%
DL_fit=DL_PF;
cumN_fit=N_PF(index_N1(index_chi_square),:);

if debug
    loglog(DL_fit,cumN_fit,'r-');hold on
    errorbar(DL,N,err_N,err_N,'rs')
    figure;semilogx(N1_PF(index_N1),chi_square,'r',N1_PF(index_N1),(min(chi_square)+2.3)*ones(size(chi_square)),'k--')
    ylim([min(chi_square)-2 min(chi_square)+10])
end

end


function [N1_fitting,N1_err,DL_fit,cumN_fit]=Chi_square_fittingLog(DL,N,err_N,N_PF,N1_PF,DL_PF,PF_slope)
if nargin<4
    N1_PF=logspace(-10,-1,1000);
    DL_PF=logspace(-3,3,500);%km
    
    for i=1:length(N1_PF)
        N_PF(i,:)=N1_PF(i)*DL_PF.^-3.2;
        PF_slope_log(i,:)=-3.2*ones(size(DL_PF));%y=ax^-b=> logy=log(a)-b*log(x)=>logy=log(a)-b*log(x)
    end
    
end
debug=0;
L=length(DL);

bin_width=diff(DL);
bin_width=[bin_width,bin_width(end)];%??????????????????

N_Dmin=interp1(DL_PF,N_PF(1,:),DL(1));
PF=interp1(DL_PF,N_PF(1,:),DL);
N_est=sum(PF.*N)/sum(PF.^2)*N1_PF(1);%derivation given in N_est.wmf in Manuscript folder.
err_N1_fit_est=(N1_PF(1)/N_Dmin)*err_N(1);


step=N1_PF(2)/N1_PF(1);
if N_est<N1_PF(1)
    fprintf(2,'N_est<min(N1_PF)\n\n');
    beep
    return;
elseif N_est>N1_PF(end)
    fprintf(2,'N_est>max(N1_PF)\n\n');
    beep
    return;
else
    index_est=log(N_est/N1_PF(1))/log(step)+1;%N_est=N1_PF(1)*step^(k-1) => k=log(N_est/N1_PF(1))/log(step)+1
    
    N_estL=N_est-50*err_N1_fit_est;
    if N_estL<N1_PF(1)
        index_L=1;
    else
        index_L=log(N_estL/N1_PF(1))/log(step)+1;
    end
    N_estU=N_est+1000*err_N1_fit_est;
    if N_estU>N1_PF(end)
        index_U=length(N1_PF);
    else
        index_U=log(N_estU/N1_PF(1))/log(step)+1;
    end
    index_N1=round(index_L:index_U);
end
if isempty(index_N1)
    index_N1=round(max(1,index_est-5):min(index_est+5,length(N1_PF)));
end
if debug
    figure;
end
err_N_log=err_N./N/log(2);
sigma=0.1;
for i=1:length(index_N1)
    PF=interp1(DL_PF,N_PF(index_N1(i),:),DL);
    slope_log=interp1(DL_PF,PF_slope_log(index_N1(i),:),DL,'next');
    weight=err_N_log.^2+slope_log.^2.*sigma.^2;%Barker 1974
    chi_square(i)=sum((log(PF)-log(N)).^2/(weight));
    
%     loglog(DL_PF,N_PF(index_N1(i),:));hold on
end


[min_chi_square,index_chi_square]=min(chi_square);
N1_fitting=N1_PF(index_N1(index_chi_square));

index=find(chi_square<min_chi_square+2.3,1,'first');
err_L=log(N1_PF(index_N1(index)));
index=find(chi_square<min_chi_square+2.3,1,'last');
err_U=log(N1_PF(index_N1(index)));
N1_err=(err_U-err_L)/2;%log
N1_err=N1_fitting*N1_err;%convert log to linear
DL_fit=DL_PF;
cumN_fit=N_PF(index_N1(index_chi_square),:);

if debug
    loglog(DL_fit,cumN_fit,'r-');hold on
    errorbar(DL,N,err_N,err_N,'rs')
    figure;plot(log(N1_PF(index_N1)),chi_square,'r',log(N1_PF(index_N1)),(min(chi_square)+2.3)*ones(size(chi_square)),'k--')
    ylim([min(chi_square)-2 min(chi_square)+10])
end

end


function [DL,cumN,err_cumN]=shp2cum(filename,Area,step,DL_pass)
%[DL,cumN,err_cumN]=shp2cumN(filename,Area,step)
if nargin<1
    close all
    filename='I:\Documents\Paper\Small craters\ArcGIS\Apollo 16 From Xiao\Craters\CRATER_Kaguya-A16.dbf';
    Area=879;%km^2
    step='pseudo-log';
end

data=shaperead(filename);
L=length(data);
D=zeros(1,L);
for i=1:L
    D(i)=data(i).Diam_km;
end

DL=sort(D,'descend');
cumN=cumsum(DL>0);

DL=DL(end:-1:1);
cumN=cumN(end:-1:1);

if isempty(DL)
   	err_cumN=[];
    cumN=[]; 
    return
end
if nargin>=3||nargin<1
    if length(step)==1
        if nargin==3
            DL_pass=1;
        end
    %     DL=1/step^round((log10(1/(min(DL)))/log10(step)));
        Num_step=floor(log10(DL_pass/min(D))/log10(step));
        DL_bin=DL_pass/step^Num_step;
        if DL_bin>DL_pass
            DL_bin=DL_pass;
        end
        i=1;
        max_D_km=max(D);
        while DL_bin<max_D_km
            i=i+1;
            DL_bin(i,1)=DL_bin(i-1,1)*step;
        end
        DL_bin(end)=[];
    else
        if strcmp(step,'pseudo-log')
            DL_bin=[];
            for i=floor(log10(min(DL))):ceil(log10(max(DL)))-1
                DL_bin=[DL_bin,[1.,1.1,1.2,1.3,1.4,1.5,1.7,2.,2.5,3.0,3.5,4.,4.5,5.,6.,7.,8.,9.] * 10^i];
            end
            index=DL_bin>min(DL);
            DL_bin=DL_bin([index(2:end) 1==1]);
        end
    end

    %------------------------直径有相同的,DL必须是升序
    [DL,index]=unique(DL);
    cumN=cumN(index);
   
    
    %------------------------
    cumN=interp1(DL,cumN,DL_bin,'next','extrap');
    DL=DL_bin;
    
    index=isnan(cumN);
    cumN(index)=[];
    DL(index)=[];
    
    if 1
        [cumN,index]=unique(fliplr(cumN));%------------------------N相同的,只要对应直径最大的
        DL=fliplr(DL);

        DL=fliplr(DL(index));
        cumN=fliplr(cumN);
    end
    
    index=cumN>0;
    cumN=cumN(index);
    DL=DL(index);
end
err_cumN=cumN.^0.5/Area;
cumN=cumN/Area;
%     for j=1:length(D1)
%         if j==1
%             N(j)=1
%         else
%             N(j)=N(j-1)+1;
%         end
%     end
    
%     N=interp1(D1,N,D2,'next','extrap');
if nargin<1
    figure;
    set(gcf,'position',[420,5,630,670])%设置figure窗口的大小
    errorbar(DL,cumN,err_cumN,'gh','MarkerSize',6);hold on
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    
    log_axis_equal
%     xlim([0.5 3]);
%     ylim([0 80])

    set(gca,'FontName','Times New Roman');
    
    label_FontSize=25;
    xlabel('Crater Diameter (km)','FontSize',label_FontSize,'FontName','Times New Roman');%设置横坐标名称和字体大小
    ylabel('Cumulative number','FontSize',label_FontSize,'FontName','Times New Roman');%设置纵坐标名称和字体大小
    set(gca,'FontSize',round(label_FontSize*0.9));%设置坐标轴刻度值的字体大小
    set(gca,'FontName','Times New Roman');
end
end
