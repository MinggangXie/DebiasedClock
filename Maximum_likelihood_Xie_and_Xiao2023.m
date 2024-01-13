function [fitresult,fitresult_referenceTarget]=Maximum_likelihood_Xie_and_Xiao2023(diam_km,Dmin_fit,Area,target_type,isDegradation,DL_fit)
%Xie and Xiao, 2023, https://doi.org/10.1016/j.epsl.2022.117963
if exist('diam_km','var')&&isstruct(diam_km)
    N_cum=diam_km.N_cum;
    err=diam_km.err;
    diam_km=diam_km.diam_km;
 
    isUnbinnedCSFD=0;
else
    isUnbinnedCSFD=1;
end
if ~exist('DL_fit','var')
    DL_fit=logspace(-3,log10(300),1000)';
end
if exist('Area','var')&&isUnbinnedCSFD==1
    N_cum=(length(diam_km):-1:1)/Area;
end
if ~exist('isDegradation','var')
    isDegradation=1;
end
if nargin<1
    dbstop if error
    close(figure(10))
    [filename,PATHNAME]=uigetfile( ...
        {'*.scc;*.shp;*.dbf;*.csfd;*.txt'},...
        'Please select data file(s) (multiple options are available)', '.\data',...
        'MultiSelect', 'on');
    if ~(ischar(filename)||iscell(filename))
        return;
    end
    if strcmp(filename(end-2:end),'txt')%has to be unbinned CSFD
        data=load([PATHNAME filename]);
        diam_km=data(:,1)';
        [D_sort,index]=sort(diam_km);
        N_cum=data(index,2)';
        Area = str2num(cell2mat(inputdlg('Area of the unit (km^{-2})')));
        dnum=(N_cum-[N_cum(2:end),0])*Area;
    elseif strcmp(filename(end-2:end),'dbf')%Robbins global dataset; has to be unbinned CSFD
        [NUM,TXT,RAW]=xlsread([PATHNAME filename]);
        diam_km=NUM(:,3)';
        [D_sort,index]=sort(diam_km);
        Area = str2num(cell2mat(inputdlg('Area of the unit (km^{2})')));
        N_cum=(length(D_sort):-1:1)/Area;
        dnum=(N_cum-[N_cum(2:end),0])*Area;
    else
        [diam_km,~,~,Area]=readdiam([PATHNAME filename]);
        D_sort=sort(diam_km);
        N_cum=(length(D_sort):-1:1)/Area;
    end
    
    target_type_list={'pure_regolith','regolith_on_megaregolith_on_rocks','regolith_on_rocks'};
    [sel,ok]=listdlg('ListString',{'pure_regolith','regolith_on_megaregolith_on_rocks','regolith_on_rocks'},...
        'Name','Target type','OKString','OK','CancelString','Cancel','SelectionMode','single','ListSize',[280 100]);
    if ok==1
        target_type=target_type_list{sel};
        %         target_type='regolith_on_rocks';%{'pure_regolith','regolith_on_rocks'};
        %     target_type='regolith_on_megaregolith_on_rocks';
        %     target_type='Serenitatis';
    else
        return;
    end
    
    figure(10);set(gcf,'position',[600 60 500 700])
    clf
    h_axes=axes;
    
    h(1)=loglog(D_sort,N_cum,'rs');hold on
    
    loglog(h_axes,DL_fit,0.01*1.54*DL_fit.^-2,'m--','linewidth',1)
    loglog(h_axes,DL_fit,0.1*1.54*DL_fit.^-2,'m--','linewidth',1)
    loglog(h_axes,DL_fit,0.06*1.54*DL_fit.^-2,'m--','linewidth',1)
    %     loglog(h_axes,DL_fit,0.04*1.54*DL_fit.^-2,'m--','linewidth',1)
    
    xlim([D_sort(1)/2,D_sort(end)*2])
    ylim([0.5*N_cum(end),2*N_cum(1)])
    
    diam_km=D_sort;
end

if ~exist('Dmin_fit','var')
    Dmin_fit = str2num(cell2mat(inputdlg('Dmin fit (km)')));
end

h_waitbar = waitbar(0,'Fitting, please wait...');
set(h_waitbar, 'WindowStyle','modal');
waitbar(0,h_waitbar);
Dsc=19;%km

N_obs_in=sum(diam_km>=Dmin_fit);
% Dmax=min(Dsc,sqrt(Area));%km; The transition from simple to complex craters is assumed to occur at 19 km
Dmax=sqrt(Area/4);%km; The transition from simple to complex craters is assumed to occur at 19 km
if Dmax<Dmin_fit
    fprintf(2,'Dmax<Dmin_fit in Maximum_likelihood_Michael2016.m\n\n')
    beep;
    return;
end
% if ~exist('DL_fit','var')
%     DL_fit=logspace(log10(min([min(diam_km)/2,1e-3])),log10(200),1000)';
index=find(DL_fit>=Dmin_fit,1,'first');
DL_fit=DL_fit/DL_fit(index)*(Dmin_fit*1.000000001);
DL_fit_EPF=DL_fit(DL_fit>=Dmin_fit&DL_fit<=Dmax);
% end

[N_Dmin_fit,dataset]=getPF(Dmin_fit,target_type);
if ~isDegradation
    dataset.N1_mpf=dataset.N1_mopf;
end
N_Dmin_fit=N_Dmin_fit*Area;%number
P_poisson=zeros(size(N_Dmin_fit));
P_log=zeros(size(N_Dmin_fit));
L=length(N_Dmin_fit);

index_used_craters=find(diam_km>=Dmin_fit);
num_used_craters=length(index_used_craters);%for unbinned craters
num_cum_max_round=round(N_cum(index_used_craters(1))*Area);
num_cum_max=N_cum(index_used_craters(1))*Area;
for i=1:L
    %     P_poisson(i)=poisspdf(N_obs_in,N_Dmin_fit(i));
    
    %         CSFD=dataset.CMPF_TDs(index,i);
    %         DSFD=dataset.DMPF_TDs(index,i);
    %         c=sum(log(DSFD./CSFD))+sum(log(bin_width))+sum(log(1:N_obs_in));
    %         if P_poisson(i)==0
    %             P_log(i)=-inf;
    %         else
    %             P_log(i)=c+log(P_poisson(i));
    %         end
    index_lessThanDmax=find(dataset.D>Dmin_fit*0.9&dataset.D<Dmax);
    if isUnbinnedCSFD
        if isDegradation
            [DMPF_TDE,~,CMPF_TDE,CMPF_TD]=get_MPF_TDE_from_MPF_TD(dataset.D(index_lessThanDmax),dataset.DMPF_TDs(index_lessThanDmax,i),'dsfd');
%             figure;
%             loglog(dataset.DL(index_lessThanDmax),dataset.CMPF_TDs(index_lessThanDmax,i),'k-');hold on
%             loglog(dataset.DL(index_lessThanDmax),CMPF_TD,'g--');hold on
%             loglog(dataset.DL(index_lessThanDmax),CMPF_TDE,'r--');
        else
            index2_R=index_lessThanDmax+1;
            DMPF_T=(dataset.CMPF_Ts(index_lessThanDmax,i)-dataset.CMPF_Ts(index2_R,i))./dataset.dD(index_lessThanDmax);
            [DMPF_TDE,~,CMPF_TDE]=get_MPF_TDE_from_MPF_TD(dataset.D(index_lessThanDmax),DMPF_T,'dsfd');
        end
        %             figure(2);loglog(dataset.D(index_lessThanDmax),dataset.DMPF_TDs(index_lessThanDmax,i),'r',dataset.DL(index_lessThanDmax),DMPF_TDE,'k--')
        %             figure(3);loglog(dataset.DL(index_lessThanDmax),dataset.CMPF_TDs(index_lessThanDmax,i),'r',dataset.DL(index_lessThanDmax),CMPF_TDE,'k--')
        C=exp(interp1(log(dataset.DL(index_lessThanDmax,1)),log(CMPF_TDE),log(Dmin_fit),'linear','extrap'));
        DSFD=exp(interp1(log(dataset.D(index_lessThanDmax,1)),log(DMPF_TDE),log(diam_km(index_used_craters)),'linear','extrap'));%DSFD and ISFD should use D rather than DL
        N_Dmin_fit(i,1)=exp(interp1(log(dataset.DL(index_lessThanDmax,1)),log(CMPF_TDE),log(Dmin_fit),'linear','extrap'));
        %             P_log(i)=N_obs_in*log(Area)+sum(log(dataset.dD(index)))-Area*C+sum(log(DSFD));
        if num_used_craters==num_cum_max_round%has to be unbinned CSFD
            n_bar=(C-CMPF_TDE(end))*Area;
            P_log(i)=-n_bar+sum(log(DSFD));
        else% may be problematic！！！
            n_bar=(C-CMPF_TDE(end))*Area;
            P_log(i)=-n_bar+sum(log(DSFD).*dnum(index_used_craters));
            %                 DSFD=DSFD*num_used_craters/num_cum_max;
        end
        
%         index=dataset.DL(:,1)>Dmin_fit;
%         C=interp1(dataset.DL(:,1),dataset.CMPF_TDs(:,i),Dmin_fit);
%         DSFD=exp(interp1(log(dataset.D(:,1)),log(dataset.DMPF_TDs(:,i)),log(diam_km(diam_km>=Dmin_fit))));
%         P_log(i)=N_obs_in*log(Area)+sum(log(dataset.dD(index)))-Area*C+sum(log(DSFD));
%         n_bar=(N_Dmin_fit(i)-interp1(dataset.DL(index,1),dataset.CMPF_TDs(index,i),Dmax));
%         
%         P_log(i)=N_obs_in*log(n_bar)-2*n_bar+sum(log(DSFD))-sum(1:N_obs_in);
    else
        if isDegradation
            [~,~,CMPF_TDE]=get_MPF_TDE_from_MPF_TD(dataset.D(index_lessThanDmax),dataset.DMPF_TDs(index_lessThanDmax,i),'dsfd');
        else
            index2_R=index_lessThanDmax+1;
            DMPF_T=(dataset.CMPF_Ts(index_lessThanDmax,i)-dataset.CMPF_Ts(index2_R,i))./dataset.dD(index_lessThanDmax);
            [~,~,CMPF_TDE]=get_MPF_TDE_from_MPF_TD(dataset.D(index_lessThanDmax),DMPF_T,'dsfd');
        end
        CSFD=exp(interp1(log(dataset.DL(index_lessThanDmax,1)),log(CMPF_TDE),log(diam_km(index_used_craters)),'linear','extrap'));%DSFD and ISFD should use D rather than DL
        chisquare(i)=sum((CSFD-N_cum(index_used_craters)).^2./err(index_used_craters).^2);
    end
    
    
    if mod(i,round(L/200))==1
        waitbar(i/L,h_waitbar)
    end
    %     P_log(i)=N_obs_in*log(Area)+sum(log(bin_width))-N_obs_in+sum(log(DSFD));
end
close(h_waitbar);
if isUnbinnedCSFD
    if 1
        %     P_log=P_log-median(P_log(P_log>-inf&~isnan(P_log)));
        P_log=P_log-median(P_log);
%         if max(P_log)>700
            P_log=P_log-max(P_log)+50;
%         end
        PDF=exp(P_log);%probability density
    else
        PDF=P_poisson;
    end

    PDF_cum=cumsum(PDF);
    PDF=PDF/PDF_cum(end);
%     PDF_cum=PDF_cum/PDF_cum(end);
    P=PDF.*[0;diff(dataset.N1_mpf)];%probability
    P_cum=cumsum(P);
    P_cum=P_cum/P_cum(end);
    [P_cum_unique,index_unique]=unique(P_cum);

    [~,index_bestFit]=max(PDF);
    N1_mpf_fit=dataset.N1_mpf(index_bestFit);

    N1_mpf_err_neg=interp1(P_cum_unique,dataset.N1_mpf(index_unique),(1-0.68)/2,'linear','extrap');
    N1_mpf_err_pos=interp1(P_cum_unique,dataset.N1_mpf(index_unique),1-(1-0.68)/2,'linear','extrap');
    N1_mopf_fit=dataset.N1_mopf(index_bestFit);
    N1_mopf_err_neg=interp1(P_cum_unique,dataset.N1_mopf(index_unique),(1-0.68)/2,'linear','extrap');
    N1_mopf_err_pos=interp1(P_cum_unique,dataset.N1_mopf(index_unique),1-(1-0.68)/2,'linear','extrap');
else
    [chisquare_min,index]=min(chisquare);

    N1_mpf_fit=dataset.N1_mpf(index);
    N1_mopf_fit=dataset.N1_mopf(index);
    idx=find(chisquare-chisquare_min<1,1,'first');%1 sigma;
    N1_mpf_err_neg=dataset.N1_mpf(idx);
    N1_mopf_err_neg=dataset.N1_mopf(idx);
    idx=find(chisquare-chisquare_min<1,1,'last');%1 sigma;
    N1_mpf_err_pos=dataset.N1_mpf(idx);
    N1_mopf_err_pos=dataset.N1_mopf(idx);
end
% index_bestFit=find(N1_mopf_fit>dataset.N1_mopf,1,'last');
CMPF_TD_best_fit=exp(interp1(log(dataset.DL),log(dataset.CMPF_TDs(:,index_bestFit)/dataset.N1_mopf(index_bestFit)*N1_mopf_fit),log(DL_fit),'linear','extrap'));
CMPF_T_best_fit=exp(interp1(log(dataset.DL),log(dataset.CMPF_Ts(:,index_bestFit)/dataset.N1_mopf(index_bestFit)*N1_mopf_fit),log(DL_fit),'linear','extrap'));



index=find(dataset.D<Dmax);
if isDegradation
    [~,~,CMPF_TDE]=get_MPF_TDE_from_MPF_TD(dataset.D(index),dataset.DMPF_TDs(index,index_bestFit)/dataset.N1_mopf(index_bestFit)*N1_mopf_fit,'dsfd');
else
    DMPF_T=(dataset.CMPF_Ts(index,index_bestFit)-dataset.CMPF_Ts(index+1,index_bestFit))./dataset.dD(index);
    [~,~,CMPF_TDE]=get_MPF_TDE_from_MPF_TD(dataset.D(index),DMPF_T/dataset.N1_mopf(index_bestFit)*N1_mopf_fit,'dsfd');
end
CEPF_best_fit=exp(interp1(log(dataset.DL(index)),log(CMPF_TDE),log(DL_fit_EPF),'linear','extrap'));

if nargin<1
    h(end+1:end+3)=loglog(DL_fit(DL_fit>=Dmin_fit),CMPF_TD_best_fit(DL_fit>=Dmin_fit),'k',DL_fit_EPF,CEPF_best_fit,'k--',DL_fit,CMPF_T_best_fit,'r-.','LineWidth',2);hold on
    [DL,cumN]=NPF(N1_mopf_fit);
    h(end+1)=loglog(DL,cumN,'g--','LineWidth',1.5);
    
    if strcmp(target_type,'pure_regolith')
        str_target='pure regolith';
    elseif strcmp(target_type,'regolith_on_megaregolith_on_rocks')
        str_target='megaregolith on bedrock';
    elseif strcmp(target_type,'regolith_on_rocks')
        str_target='regolith on bedrock';
    end
    str_leg={'measured SFD',['CMPF_{TD} (' str_target ')'],'Bestfit CMPF_{TDE}',sprintf('CMPF_{T} = %.2g_{-%.2g}^{+%.2g} km^{-2}',N1_mpf_fit,N1_mpf_fit-N1_mpf_err_neg,N1_mpf_err_pos-N1_mpf_fit),'NPF'};
    title(filename(1:end-4));
    
    FontSize=15;
    xlabel('Crater diameter (km)','FontSize',FontSize,'FontName','Times New Roman')
    ylabel('Cumulative crater density (km^{-2})','FontSize',FontSize,'FontName','Times New Roman')
    %--------------------------------------inset-------------------------------------------------
    h_inset=axes('Position',[0.15,0.25,0.35,0.25]);
    index=PDF>0.01*max(PDF);
    plot(h_inset,dataset.N1_mpf(index),PDF(index));
    set(h_inset,'Xscale','log','box','off','ytick',[],'ycolor','w','FontSize',FontSize*0.5);
    xlabel('Crater density (km^{-2})','FontSize',FontSize*0.7,'FontName','Times New Roman')
    %---------------------------------------------------------------------------------------
    %     plot(N_Dmin_fit(index),PDF(index))
end
if isUnbinnedCSFD
    fitresult.PDF=PDF;
else
    fitresult.chisquare=chisquare;
end
fitresult.N_Dmin_fit=N_Dmin_fit;
fitresult.N1_mpf_fit=N1_mpf_fit;
fitresult.N1_mpf_fit=N1_mpf_fit;
fitresult.N1_mpf_dataset=dataset.N1_mpf;
fitresult.N1_mopf_dataset=dataset.N1_mopf;

fitresult.N1_mpf_err_neg=N1_mpf_err_neg;
fitresult.N1_mpf_err_pos=N1_mpf_err_pos;
fitresult.N1_mopf_fit=N1_mopf_fit;
fitresult.N1_mopf_err_neg=N1_mopf_err_neg;
fitresult.N1_mopf_err_pos=N1_mopf_err_pos;
fitresult.DL_fit=DL_fit;
fitresult.D_fit=DL_fit*(DL_fit(2)/DL_fit(1))^0.5;
fitresult.CMPF_TD_best_fit=CMPF_TD_best_fit;
fitresult.CMPF_T_best_fit=CMPF_T_best_fit;
fitresult.CEPF_best_fit=CEPF_best_fit;
fitresult.DL_fit_EPF=DL_fit_EPF;
fitresult.D_fit_EPF=DL_fit_EPF*(DL_fit_EPF(2)/DL_fit_EPF(1))^0.5;

if nargout==2||nargin<1
    target_type='reference_target';
    
    [~,dataset]=getPF(Dmin_fit,target_type);
    if isUnbinnedCSFD
        [~,index_bestFit]=max(PDF);
        N1_mpf_fit=dataset.N1_mpf(index_bestFit);

        N1_mpf_err_neg=interp1(P_cum_unique,dataset.N1_mpf(index_unique),(1-0.68)/2,'linear','extrap');
        N1_mpf_err_pos=interp1(P_cum_unique,dataset.N1_mpf(index_unique),1-(1-0.68)/2,'linear','extrap');
        N1_mopf_fit=dataset.N1_mpf(index_bestFit);
        N1_mopf_err_neg=interp1(P_cum_unique,dataset.N1_mpf(index_unique),(1-0.68)/2,'linear','extrap');
        N1_mopf_err_pos=interp1(P_cum_unique,dataset.N1_mpf(index_unique),1-(1-0.68)/2,'linear','extrap');
    else
        [chisquare_min,index]=min(chisquare);

        N1_mpf_fit=dataset.N1_mpf(index);
        idx=find(chisquare-chisquare_min<1,1,'first');%1 sigma;
        N1_mpf_err_neg=dataset.N1_mpf(idx);
        idx=find(chisquare-chisquare_min<1,1,'last');%1 sigma;
        N1_mpf_err_pos=dataset.N1_mpf(idx);
    end
%     index_bestFit=find(N1_mopf_fit>dataset.N1_mpf,1,'last');
    %     CMPF_TD_best_fit=exp(interp1(log(dataset.DL),log(dataset.CMPF_Ts(:,index_bestFit)),log(DL_fit)));
    CMPF_T_best_fit=exp(interp1(log(dataset.DL),log(dataset.CMPF_Ts(:,index_bestFit)),log(DL_fit),'linear','extrap'));
    if nargin<1
        h(end+1)=loglog(h_axes,DL_fit,CMPF_T_best_fit,'b:','LineWidth',2);hold on
        
        age=debiased_density2age(N1_mpf_fit);
        str_age=sprintf('%.2g',round(age,3,'significant'));
        str_leg={str_leg{:},sprintf('Reference CMPF_{T} with debiased N(1) = %.2g_{-%.2g}^{+%.2g} km^{-2} \n and AMA = %.3g_{-%.2g}^{+%.2g} Ga',N1_mpf_fit,N1_mpf_fit-N1_mpf_err_neg,N1_mpf_err_pos-N1_mpf_fit,...
        debiased_density2age(N1_mpf_fit),debiased_density2age(N1_mpf_fit)-debiased_density2age(N1_mpf_err_neg),debiased_density2age(N1_mpf_err_pos)-debiased_density2age(N1_mpf_fit))};
        legend(h,str_leg);
        
        
        
        xlim(h_axes,[min(diam_km)*0.8,max(max(diam_km),1)*1.2])
        %         ylim([1e-7,2e6])
    end
    fitresult_referenceTarget.N1_mpf_fit=N1_mpf_fit;
    
    fitresult_referenceTarget.N1_mpf_err_neg=N1_mpf_err_neg;
    fitresult_referenceTarget.N1_mpf_err_pos=N1_mpf_err_pos;
    fitresult_referenceTarget.N1_mopf_fit=N1_mopf_fit;
    fitresult_referenceTarget.N1_mopf_err_neg=N1_mopf_err_neg;
    fitresult_referenceTarget.N1_mopf_err_pos=N1_mopf_err_pos;
    fitresult_referenceTarget.DL_fit=DL_fit;
    %     fitresult_referenceTarget.CMPF_TD_best_fit=CMPF_TD_best_fit;
    fitresult_referenceTarget.CMPF_T_best_fit=CMPF_T_best_fit;
    fitresult_referenceTarget.N1_mpf_dataset=dataset.N1_mpf;
end


%--------------------------------------------------------------------------------------------------------------------------------------
function [N_Dmin_fit,dataset]=getPF(Dmin_fit,target_type,PF_name)%D,km;
if nargin<3
    PF_name='XPF';
end
switch PF_name
    case 'XPF'
        %         if strcmp(target_type,'Serenitatis')
        %             dataset=load('..\MOPF\version2\data\MOPF_regolith_on_fractured_rocks.mat');
        %             dataset.CMOPFs=dataset.CMOPFs*100;
        %             dataset.CMPFs=dataset.CMPFs*100;
        %             dataset.DMPF_TDs=dataset.DMPF_TDs*100;
        %         else
        filename=['..\MPF\version3\data\MPF_TD_', target_type '.mat'];
        dataset=load(filename);
        fprintf(2,'loading data from version3 \n\n')
        %         end
    case 'NPF'
        published_time=2001;
        [DL,cumN]=NPF(N1,published_time);
        dataset.DL=DL;
        dataset.CMOPFs=cumN;
        
end

step=dataset.DL(2)/dataset.DL(1);
L_D=length(dataset.DL);
N_Dmin_fit=zeros(length(dataset.CMPF_Ts(1,:)),1);
idx=round(log(Dmin_fit/dataset.DL(1))/log(step)+1);%a(1)*step_a^(n-1)=x;
index=max(1,idx-1):min(idx+1,L_D);

for i=1:length(dataset.CMPF_Ts(1,:))
    if strcmp(target_type,'reference_target')
        N_Dmin_fit(i,1)=exp(interp1(log(dataset.DL(index,1)),log(dataset.CMPF_Ts(index,i)),log(Dmin_fit),'linear','extrap'));
    else
        N_Dmin_fit(i,1)=exp(interp1(log(dataset.DL(index,1)),log(dataset.CMPF_TDs(index,i)),log(Dmin_fit),'linear','extrap'));
    end
%     figure;
%     loglog(dataset.DL,dataset.CMPF_TDs(:,i),'k-');hold on
%     loglog(dataset.DL,CMPF_TD,'g--');hold on
%     loglog(dataset.DL,CMPF_TDE,'r--');
end

function [DL,cumN]=NPF(N1,published_time)
%[DL,cumN]=NPF(t,Drange_km,step,published_time)
%The unit of t is Ga
%step=DR/DL
%----------------
if nargin<1
    Drange_km=[0.01,100];
    t=3;
    step=2^(1/4);
    step=1/0.8;
end
if nargin<4
    published_time=2001;
end
DL=logspace(-2,log10(3e2),1000);
NPF1=zeros(size(DL));
if published_time==2001
    aj=[-3.0876 -3.557528 0.781027 1.021521 -0.156012 -0.444058 0.019977 0.086850 -0.005874 -0.006809 8.25e-4 5.54e-5];
    
    %     N1=5.44e-14*(exp(6.93*t)-1)+8.38e-4*t;
    aj(1)=log10(N1);
    for i=1:length(aj)
        NPF1=NPF1+aj(i).*(log10(DL)).^(i-1);
    end
elseif published_time==1983
    aj=[-3.0768 % 0
        -3.6269 % 1
        0.4366 % 2
        0.7935 % 3
        0.0865 % 4
        -0.2649 % 5
        -0.0664 % 6
        0.0379 % 7
        0.0106 % 8
        -0.0022 % 9
        -5.18e-4 % 10
        3.97e-5];% 11
    
    %     N1=5.44e-14*(exp(6.93*t)-1)+8.38e-4*t;
    aj(1)=log10(N1);
    for i=1:length(aj)
        NPF1=NPF1+aj(i).*(log10(DL)).^(i-1);
    end
end
cumN=10.^NPF1;
if nargin<1
    loglog(DL,cumN,'r-');
end


function [DMPF_TDE,IMPF_TDE,CMPF_TDE,CMPF_TD]=get_MPF_TDE_from_MPF_TD2(edf_D,DMPF_TD,type)
if nargin<3
    type='dsfd';
end
sigma_var=0.1;
step=edf_D(2)/edf_D(1);
edf_DR=edf_D*step^0.5;
edf_DL=edf_D*step^-0.5;

bin_width=edf_DR-edf_DL;

DMPF_TDE=zeros(size(edf_D));
if strcmp(type, 'dsfd')
    MOPF_isfd=DMPF_TD.*bin_width;
    D=edf_D;
elseif strcmp(type, 'diameters')
    D=DMPF_TD;
    MOPF_isfd=ones(size(D));
end
sigma=sigma_var*D;

for i=1:length(MOPF_isfd)
    index=find(D(i)-4*sigma(i)<edf_D*step&edf_D<=D(i)+4*sigma(i));
    y=1/(sqrt(2*pi)*sigma(i))*exp(-0.5*((edf_D(index)-D(i))/sigma(i)).^2);
    %     y=y/sum(y.*bin_width(index))*MOPF_isfd(i);%
    y=y*MOPF_isfd(i);%不能通过除以sum(y.*bin_width(index))来归一化，因为有可能某些情况是在横坐标范围外的，比如L=0.1,而一个0.1直径的撞击坑，有一半概率是小于L的
    DMPF_TDE(index)=DMPF_TDE(index)+y;
    %     p(i,index)=y;
    %     loglog(edf.DL(index),y,'r--');hold on
end
IMPF_TDE=DMPF_TDE.*bin_width;
CMPF_TDE=cumsum(IMPF_TDE(end:-1:1));
CMPF_TDE(end:-1:1)=CMPF_TDE;

CMPF_TD=cumsum(MOPF_isfd,1,'reverse');

function [DMPF_TDE,IMPF_TDE,CMPF_TDE,CMPF_TD]=get_MPF_TDE_from_MPF_TD(edf_D,DMPF_TD,type)
if nargin<1
    close all
    edf_D=logspace(-1,1,100);
    DMPF_TD=edf_D.^-4;
end
if nargin<3
    type='dsfd';
end
sigma_var=0.03;
step=edf_D(2)/edf_D(1);
edf_DR=edf_D*step^0.5;
edf_DL=edf_D*step^-0.5;

bin_width=edf_DR-edf_DL;

if strcmp(type, 'dsfd')
    MOPF_isfd=DMPF_TD.*bin_width;
    D=edf_D;
elseif strcmp(type, 'diameters')
    disp('undefined ''type of diameters''in get_MPF_TDE_from_MPF_TD')
    beep
    return;
    D=DMPF_TD;
    MOPF_isfd=ones(size(D));
end

n_bin=round(log(max(edf_DL)/min(edf_DL))/log(1+sigma_var/3));
DL_fine=logspace(log10(min(edf_DL)),log10(max(edf_DL)),n_bin);
D_fine=DL_fine*(DL_fine(2)/DL_fine(1))^0.5;
DR_fine=DL_fine*DL_fine(2)/DL_fine(1);
DMPF_TD_fine=exp(interp1(log(edf_D),log(DMPF_TD),log(D_fine),'linear','extrap'));
IMPF_fine=DMPF_TD_fine.*(DR_fine-DL_fine);
DMPF_TDE_fine=zeros(size(IMPF_fine));

sigma=sigma_var*D_fine;
for i=1:n_bin
    index=find(D_fine(i)-4*sigma(i)<DR_fine&DL_fine<=D_fine(i)+4*sigma(i));
    y=1/(sqrt(2*pi)*sigma(i))*exp(-0.5*((D_fine(index)-D_fine(i))/sigma(i)).^2);
    %     y=y/sum(y.*bin_width(index))*MOPF_isfd(i);%
    y=y*IMPF_fine(i);%不能通过除以sum(y.*bin_width(index))来归一化，因为有可能某些情况是在横坐标范围外的，比如L=0.1,而一个0.1直径的撞击坑，有一半概率是小于L的
    DMPF_TDE_fine(index)=DMPF_TDE_fine(index)+y;
    %     p(i,index)=y;
    %     loglog(edf.DL(index),y,'r--');hold on
end
if length(MOPF_isfd(1,:))>length(MOPF_isfd(:,1))
    CMPF_TD=cumsum(MOPF_isfd,2,'reverse');
else
    CMPF_TD=cumsum(MOPF_isfd,1,'reverse');
end
DMPF_TDE=exp(interp1(log(D_fine),log(DMPF_TDE_fine),log(D),'linear','extrap'));
IMPF_TDE=DMPF_TDE.*bin_width;
% CMPF_TDE=cumsum(IMPF_TDE(end:-1:1));
% CMPF_TDE(end:-1:1)=CMPF_TDE;

if length(IMPF_TDE(1,:))>length(IMPF_TDE(:,1))
    CMPF_TDE=cumsum(IMPF_TDE,2,'reverse');
else
    CMPF_TDE=cumsum(IMPF_TDE,1,'reverse');
end
if nargin<1
    figure;
    subplot(2,1,1)
    loglog(edf_D,CMPF_TD,'k-');hold on
    loglog(edf_D,CMPF_TDE,'r--');
    subplot(2,1,2)
    loglog(edf_D,CMPF_TDE./CMPF_TD,'r-');
end



