function handles = main(handles,isjustPlot)
if nargin<2
    isjustPlot=0;
end
index_selected= get(handles.listbox1,'value'); %index的值代表我们选的是第几个选项
if isempty(index_selected)
    return;
end
filename=get(handles.listbox1,'string'); %获取列表中的所有数据
if isempty(filename)
    return;
end
if ischar(filename)
    filename={filename};
end
Ylim=[inf -inf];
% axes(handles.axes1);
% hold off
cla %reset

data=handles.data;
data.index_selected=index_selected;
pathname=data.pathname;

PF_type='NPF';
DL_min=inf;
DL_max=0;
for i=1:length(index_selected)
    temp=filename{index_selected(i)};
    filefullname=[pathname{index_selected(i)} filename{index_selected(i)}];
    if ~exist(filefullname,'file')
        msgbox('File don''t exist');
        return;
    end
    if strcmp(temp(end-3:end),'csfd')
        temp=load(filefullname);
        DL=temp(:,1);
        cumN=temp(:,2);
        err_cumN=temp(:,3);
        Area=temp(1,4);
    else
        [DL,cumN,err_cumN,Area]=readdiam(filefullname);
    end
    if DL(1)<DL_min
        DL_min=DL(1);
    end
    if DL(end)>DL_max
        DL_max=DL(end);
    end


    SFD_method=data.SFD_method{index_selected(i)};
    
    hold on;
    switch SFD_method
        case 'Unbinned CSFD'
            %------------------------直径有相同的,DL必须是升序
            h_sfddata(i)=errorbar(DL,cumN,err_cumN,err_cumN,[data.MarkerColor{index_selected(i)},data.MarkerType{index_selected(i)}],'MarkerSize',data.MarkerSize{index_selected(i)});
            set(gca,'Yscale','log','Xscale','log')
        case 'Xie 2021 (bining_CSFD)'
            %------------------------直径有相同的,DL必须是升序
            [DL2,index]=unique(DL);
            cumN2=cumN(index);
            %------------------------
            DL3=logspace(log10(min(DL)),log10(max(DL)),20)';
            cumN=interp1(DL2,cumN2,DL3,'next');
            err_cumN=(cumN*Area).^0.5/Area;
            h_sfddata(i)=errorbar(DL3,cumN,(cumN/Area).^0.5,(cumN/Area).^0.5,0.1*DL3,0.1*DL3,[data.MarkerColor{index_selected(i)},data.MarkerType{index_selected(i)}],'MarkerSize',data.MarkerSize{index_selected(i)});hold on
            set(gca,'Yscale','log','Xscale','log')
        case 'Robbins et al. (2018) (Default)'
            file_craters=dir(filefullname);
            file_craters.datenum;
            file_mat_name=['temp\' temp(1:end-4) '.mat'];
            file_mat=dir(file_mat_name);
            if isempty(file_mat)||file_mat.datenum<file_craters.datenum
                edf=SFD_robbins2018(DL,Area);
                save(file_mat_name, '-struct', 'edf');
                %                 edf=uncertainty(diam,Area,edf_diam');
            else
                edf=load(file_mat_name);
            end
            loglog(edf.diam, edf.erro_csfd_1s_pos/Area, 'g--','LineWidth', .5 ); hold on% + 1s error
            loglog( edf.diam, edf.erro_csfd_1s_neg/Area, 'g--','LineWidth', .5 ); % - 1s error
            loglog( edf.diam, edf.csfd/Area,'g','MarkerSize',8,'LineWidth',2);    % data
            %             uncertainty(DL,Area);
        case 'pseudo-log (Neukum, 1983)'%'pseudo-log'
            
            DL_bin=[];
            for j=floor(log10(min(DL))):ceil(log10(max(DL)))-1
                DL_bin=[DL_bin,[1.,1.1,1.2,1.3,1.4,1.5,1.7,2.,2.5,3.0,3.5,4.,4.5,5.,6.,7.,8.,9.] * 10^j];
            end
            index=DL_bin>min(DL);
            DL_bin=DL_bin([index(2:end) 1==1]);
            
            %------------------------直径有相同的,DL必须是升序
            [DL2,index]=unique(DL);
            cumN=cumN(index);
            %------------------------
            cumN=interp1(DL2,cumN,DL_bin,'next','extrap');
            DL2=DL_bin;
            
            index=isnan(cumN);
            cumN(index)=[];
            DL2(index)=[];
            
            if 1
                [cumN,index]=unique(fliplr(cumN));%------------------------N相同的,只要对应直径最大的
                DL2=fliplr(DL2);
                
                DL2=fliplr(DL2(index));
                cumN=fliplr(cumN);
            end
            
            index=cumN>0;
            cumN=cumN(index);
            DL2=DL2(index);
            err_cumN=(cumN*Area).^0.5/Area;
            %-----------------Plotting Marker--------------------
            warning off
            h_sfddata(i)=errorbar(DL2,cumN,err_cumN,[data.MarkerColor{index_selected(i)},data.MarkerType{index_selected(i)}],'MarkerSize',data.MarkerSize{index_selected(i)});hold on
%             h_sfddata(i)=errorbar(DL2,cumN,(cumN/Area).^0.5,(cumN/Area).^0.5,0.1*DL2*(DL2(2)/DL2(1))^0.5,0.1*DL2*(DL2(2)/DL2(1))^0.5,[data.MarkerColor{index_selected(i)},data.MarkerType{index_selected(i)}],'MarkerSize',data.MarkerSize{index_selected(i)});hold on
            drawnow% "waring off" does not work if errorbarloglog is not completed
            warning on
            %-----------------Show Marker settings--------------------
            set(gca,'Yscale','log','Xscale','log')
    end
    %---------------------Show setting from handles.data----------------------------------------------------------
    String=get(handles.target_type,'String');
    for j=1:length(String)
        if strcmp(data.target_type{index_selected(i)},String{j})
            set(handles.target_type,'value',j);%不能使用set(handles.target_type,'String',String{j})，否则会改变String的值，并且从cell变为char
        end
    end
    
    String=get(handles.Marker,'String');
    for j=1:length(String)
        temp=String{j};
        if strcmp(data.MarkerType{index_selected(i)},temp(1))
            set(handles.Marker,'value',j);
        end
    end
    
    String=get(handles.MarkerSize,'String');
    for j=1:length(String)
        if data.MarkerSize{index_selected(i)}==str2double(String)
            set(handles.MarkerSize,'value',j);
        end
    end
    
    String=get(handles.MarkerColor,'String');
    for j=1:length(String)
        temp=String{j};
        if strcmp(data.MarkerColor{index_selected(i)},temp(1))
            set(handles.MarkerColor,'value',j);
        end
    end
    data.datasymbol_handle{index_selected(i)}=h_sfddata(i);
    if isjustPlot
        if i<length(index_selected)
            continue
        else
            handles.data=data;
            return;
        end
    end
    %-----------------fitting--------------------
    if min(cumN)<Ylim(1)
        Ylim(1)=min(cumN);
    end
    if max(cumN)>Ylim(2)
        Ylim(2)=max(cumN);
    end
    
    
    switch handles.data.Fitting_method{index_selected(i)}
        case 'Xie and Xiao (2023)'
%             [N1_fitting,P,N1_err_neg,N1_err_pos,DL_fit,PF_best_fit]=Maximum_likelihood_Michael2016(DL,data.Dmin{index_selected(i)},Area);
            target_type=data.target_type{index_selected(i)};
            [fitresult,fitresult_referenceTarget]=Maximum_likelihood_Xie2023(DL,data.Dmin{index_selected(i)},Area,target_type);
            index=fitresult.DL_fit_EPF<inf;
            h_fit(i)=loglog(fitresult.DL_fit_EPF(index),fitresult.CEPF_best_fit(index),[data.LineColor{index_selected(i)},data.LineType{index_selected(i)}],'LineWidth',data.LineWidth{index_selected(i)});hold on
            index=fitresult_referenceTarget.DL_fit>=0;
            loglog(fitresult_referenceTarget.DL_fit(index),fitresult_referenceTarget.CMPF_T_best_fit(index),[data.LineColor{index_selected(i)},'-.'],'LineWidth',data.LineWidth{index_selected(i)});hold on
            
            N1_fitting=fitresult_referenceTarget.N1_mpf_fit;
            N1_err_pos=fitresult_referenceTarget.N1_mpf_err_pos;
            N1_err_neg=fitresult_referenceTarget.N1_mpf_err_neg;
%             [N1_fitting,N1_err,DL_fit,PF_best_fit]=Maximum_likelihood_estimation(DL(DL>data.Dmin{index_selected(i)}),Area);
        case 'Poisson (Michael)'
%             [N1_fitting,P,N1_err_neg,N1_err_pos,DL_fit,PF_best_fit]=Maximum_likelihood_Michael2016(DL,data.Dmin{index_selected(i)},Area);
%             target_type=data.target_type{index_selected(i)};
%             [N1_fitting,N1_err,DL_fit,PF_best_fit]=Maximum_likelihood_estimation(DL(DL>data.Dmin{index_selected(i)}),Area);
        case 'Cumulative fit'%-----------------Cumulativefit--------------------
            index_fit=find(DL>=data.Dmin{index_selected(i)}&DL<=data.Dmax{index_selected(i)});
            if isempty(index_fit)
                msgbox('Diameter range used for fitting is inappropriate!')
                return;
            else
                
                %----------------------PF--------------------------------------
                
                PF_type=get(handles.PF,'String');
                index_PF=get(handles.PF,'value');
                PF_type=PF_type{index_PF};%
                
                
                
                [N1_fitting,N1_err,DL_fit,PF_best_fit]=csfd2density(DL(index_fit),cumN(index_fit),err_cumN(index_fit),Area,PF_type);
                
                %             edf=SFD_robbins2018(DL,Area);
                
            end
        case 'None'
            %         str_leg{i}=sprintf('%s, area=%.3fkm^2',data.legendname{index_selected(i)},Area);
            str_N1=sprintf('%0.2e',interp1(DL,cumN,1,'PCHIP','extrap'));
            str_N1_err=sprintf('%0.2e',sqrt(interp1(DL,cumN,1,'PCHIP','extrap')/Area));
            str_leg{i}=sprintf('%s, N(1)=%s±%skm^{-2}',data.legendname{index_selected(i)},e2ten(str_N1),e2ten(str_N1_err));
    end
    
    %----------------------ChronologyFunction Function--------------------------------------
    ChronologyFunction=get(handles.ChronologyFunction,'String');
    index_ChronologyFunction=get(handles.ChronologyFunction,'value');
    ChronologyFunction=ChronologyFunction{index_ChronologyFunction};

    systematic_error_probability_density_function='Gaussian distribution';
    b=-(log(fitresult.CEPF_best_fit(2))-log(fitresult.CEPF_best_fit(1)))/(log(fitresult.D_fit_EPF(2))-log(fitresult.D_fit_EPF(1)));
    beta=data.systematic_error.beta{index_selected(i)};
    lambda=Area*exp(interp1(log(fitresult.DL_fit_EPF),log(fitresult.CEPF_best_fit),log(data.Dmin{index_selected(i)}),'linear','extrap'));
    [x,pdf]=pdf_of_sum(lambda,beta,b,systematic_error_probability_density_function);
    confidence_interval=0.9545;%1 sigma, 0.6827;2 sigma, 0.9545;1 sigma, 0.9973;
    idx=interp1([0;fitresult.N_Dmin_fit],0:length(fitresult.N_Dmin_fit),x/Area,'linear','extrap');
    x_referenceTarget=interp1(0:length(fitresult_referenceTarget.N1_mpf_dataset),[0;fitresult_referenceTarget.N1_mpf_dataset],idx,'linear','extrap')*Area;
    [lower_bound_1sigma,upper_bound_1sigma]=get_uncertainty(confidence_interval,x_referenceTarget,pdf);
    N1_err_pos=upper_bound_1sigma/Area;
    N1_err_neg=lower_bound_1sigma/Area;
    
    age(1)=density2age(N1_fitting,ChronologyFunction);
    age(2)=density2age(N1_err_pos,ChronologyFunction);
    age(3)=density2age(N1_err_neg,ChronologyFunction);

    %-----------------Plotting Line--------------------
    if ~strcmp(handles.data.Fitting_method{index_selected(i)},'Xie and Xiao (2023)')
        h_fit(i)=loglog(DL_fit,PF_best_fit,[data.LineColor{index_selected(i)},data.LineType{index_selected(i)}],'LineWidth',data.LineWidth{index_selected(i)});
    end
    %-----------------Show settings--------------------
    String=get(handles.line,'String');
    for j=2:length(String)
        temp=String{j};
        if strcmp(data.LineType{index_selected(i)},temp(1:2))
            set(handles.line,'value',j);
        end
    end

    String=get(handles.LineWidth,'String');
    for j=2:length(String)
        temp=String{j};
        if data.LineWidth{index_selected(i)}==str2double(temp(1:3))
            set(handles.LineWidth,'value',j);
        end
    end

    String=get(handles.LineColor,'String');
    for j=2:length(String)
        temp=String{j};
        if strcmp(data.LineColor{index_selected(i)},temp(1))
            set(handles.LineColor,'value',j);
        end
    end
    %--------------------------Legend------------------
    age=roundsd([age(1),age(2)-age(1),age(1)-age(3)],2);
    N1_fitting=roundsd(N1_fitting,2);

    if age(1)>=1
        age_unit='Ga';
    elseif age(1)>=1e-3
        age=age*1000;
        age_unit='Ma';
    elseif age(1)>=1e-6
        age=age*1e6;
        age_unit='Ka';
    elseif age(1)>=1e-9
        age=age*1e9;
        age_unit='a';
    end

    if get(handles.N1,'value')
        str_N1=sprintf('debiased N(1)={%s} _{-%s}^{+%s} km^{-2}, ',formatfloat(N1_fitting),formatfloat(N1_fitting-N1_err_neg),formatfloat(N1_err_pos-N1_fitting));
    else
        str_N1=[];
    end
    if age(1)>=10
        str_leg{i}=sprintf('  %s, %sAMA=%.0f_{-%.1f}^{+%.1f}%s, 2\\sigma, %g craters',data.legendname{index_selected(i)},str_N1,age(1),age(3),age(2),age_unit,length(DL(DL>data.Dmin{index_selected(i)})));
    else
        str_leg{i}=sprintf('  %s, %sAMA=%.1f_{-%.2f}^{+%.2f}%s, 2\\sigma, %g craters',data.legendname{index_selected(i)},str_N1,age(1),age(3),age(2),age_unit,length(DL(DL>data.Dmin{index_selected(i)})));
    end
    data.legend.text{index_selected(i)}=str_leg{i};
    data.bestfit_handle{index_selected(i)}=h_fit(i);
end
c1=10^(0.05*log10(Ylim(2)/Ylim(1)));
c2=10^(0.3*log10(Ylim(2)/Ylim(1)));
Ylim=[Ylim(1)/c1 Ylim(2)*c2];
ylim(Ylim);
xlim([DL_min*0.75 DL_max*1.5])

%----------------------Equilibrium function----------------------
Xlim=get(gca,'xlim');
if get(handles.Equilibrium_checkbox,'value')
    Deq=logspace(-10,5,1000);
    loglog(Deq,0.15*Deq.^-1.83,'k:','LineWidth',0.5);
    loglog(Deq,0.0924*Deq.^-1.83,'k:','LineWidth',0.5);
    loglog(Deq,0.015*Deq.^-1.83,'k:','LineWidth',0.5);
end
% text(10^mean(log10(Xlim)),10^mean(log10(Ylim)),'Equilibrium distribution','Rotation',-55)
%----------------------legend----------------------
% for i=1:length(str_leg)
%     str_leg{i}=['  ' str_leg{i}];
% end
drawnow;pause(0.01);
[h,object_h,OUTH,OUTM]=legend(handles.axes1,h_fit,str_leg);%[lgd,icons,plots,txt] = legend(___) additionally returns the objects used to create the legend icons, the objects plotted in the graph, and an array of the label text. Note: This syntax is not recommended. It creates a legend that does not support all graphics features. Instead, use the lgd = legend(__) syntax to return the legend object and set Legend Properties. (Emphasis mine)
% h=legend(handles.axes1,h_fit,str_leg);
% set(h,'FontUnits','normalized');
% warning off
% p=get(h,'position');
% p(1)=p(1)+0.042;
% set(h,'position',p);

% set(h,'Fontsize',11);
set(h,'box','off');

% objhdl = findobj(gcf,'Type','line');
objhdl= findobj(object_h, 'type', 'line');  %Find objects of legend of type line.
h_text= findobj(object_h, 'type', 'text');  %Find objects of legend of type text.
for i=1:length(objhdl)/2
    set(objhdl((i-1)*2+1),'Marker',data.MarkerType{index_selected(i)},'MarkerEdgeColor',data.MarkerColor{index_selected(i)},'LineWidth',1);
    %     set(h_text(i),'Fontsize',11);
    %     set(OUTH(i),'Marker',data.MarkerType{index_selected(i)},'MarkerEdgeColor',data.MarkerColor{index_selected(i)});
    %     set(h_fit(i),'color',data.LineColor{index_selected(i)},'LineStyle',data.LineType{index_selected(i)},'LineWidth',data.LineWidth{index_selected(i)},'Marker', 'none')
end

data.legend.handle=h_fit;

if length(index_selected)==1
    set(handles.Dmin,'string',num2str(data.Dmin{index_selected(1)}));
    set(handles.Dmax,'string',num2str(data.Dmax{index_selected(1)}));
else
    set(handles.Dmin,'string','N/A');
    set(handles.Dmax,'string','N/A');
end
handles.data=data;

