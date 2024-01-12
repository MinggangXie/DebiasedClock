function [DL,cumN,err_cumN,Area]=readdiam(filename)
%[DL,cumN,err_cumN,Area]=readdiam(filename)
%read data (.diam,.scc) exported from crater tool
%cumN is the cumulative number of craters per km^2
%err_cumN is the error of cumN
%filename: .diam
%Area is the area used for crater counting
if nargin<1
    filename= 'data\Xiao\Copernicus-Area 1-4.diam';
    filename= 'data\Xiao\Orientale.scc';
    Area=1;
end
isdiamfile=filename(end)=='m';
    fid=fopen(filename,'r');
    if isdiamfile
        while 1
            temp=fgetl(fid);
            if length(temp)>=10
                if strcmp(temp(1:4),'area')
                    Area=str2num(temp(8:end));
                elseif strcmp(temp(1:4),'Area')
                    Area=str2num(temp(15:end));
                end
                if strcmp(temp(1:10),'crater = {')||strcmp(temp(1:10),'#diameter,')
                    break;
                end
            end
        end
    else
        while 1
            temp=fgetl(fid);
            if length(temp)>=10
                if strcmp(temp(1:10),'Total_area')
                    index1=find(temp=='=');
                    index2=find(temp=='<');
                    Area=str2num(temp(index1+1:index2-1));
                end
                if strcmp(temp(1:10),'crater = {')
                    break;
                end
            end
        end
    end

    temp=fgetl(fid);
    j=1;
    while 1
        if isempty(temp)
            continue;
        elseif temp(1)=='}'
            break;
        end
        index=strfind(temp,char(9));
        num=str2num(temp(1:index(1)-1));
        
        data(j,1)=num(1);
        
        j=j+1;
        temp=fgetl(fid);
    end
    fclose(fid);
    
    D=data(:,1);
    DL=sort(D,'descend');
    cumN=cumsum(DL>0);
    
    DL=DL(end:-1:1);
    cumN=cumN(end:-1:1);
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
    loglog(DL,cumN,'ro','MarkerSize',10);hold on

%     xlim([0.5 3]);
%     ylim([0 80])

    set(gca,'FontName','Times New Roman');

    label_FontSize=25;
    xlabel('Crater Diameter (km)','FontSize',label_FontSize,'FontName','Times New Roman');%设置横坐标名称和字体大小
    ylabel('Cumulative number','FontSize',label_FontSize,'FontName','Times New Roman');%设置纵坐标名称和字体大小
    set(gca,'FontSize',round(label_FontSize*0.9));%设置坐标轴刻度值的字体大小
    set(gca,'FontName','Times New Roman');
end
