function [systematic_error_diameter,precision_diameter]=calibration(filenames,isRobbins_2014_ESPL_file)

data_pre=[];
N=length(filenames);
for i=1:N
    if exist('isRobbins_2014_ESPL_file','var')&&i==2
        data_current=load(filenames{i});
        data_current=data_current(:,[3 2 1]);
    else
        data_current=xlsread(filenames{i});
    end
    %         [DL,cumN,err_cumN,Area]=readdiam(filename)
    N_data_current=length(data_current(:,1));
    if i==1
        [diameters,idx]=sort(data_current(:,1));
        data_current=data_current(idx,:);
        data_new=data_current;
    elseif (length(data_current(:,1))<length(data_pre(:,1)))
        data_current(end+1:length(data_pre(:,1)),:)=0;
        data_new=zeros(size(data_current));
    end
    %     coordinates_lon_lat{i}=data_current(:,2:3);
    for j=1:i-1
        for k=1:N_data_current
            index=find(diameters(:,j)>data_current(k,1)*0.6&diameters(:,j)<data_current(k,1)*1.5);
            
            coordinates_lon_lat=data_pre(:,(j-1)*3+(2:3));
            d=distance(data_current(k,3),data_current(k,2),coordinates_lon_lat(index,2),coordinates_lon_lat(index,1),1737.4);
            
            min(d)
            
            idx=find(d<data_current(k,1)*0.5);
            FoundCrater=d(idx);
            
            if isempty(FoundCrater)
                data_new(end+1,:)=data_current(k,:);
                %                 data_current(k,:)=0;
                %                 diameters(k,i)=0;
                diameters(end+1,i)=data_current(k,1);
            elseif length(FoundCrater)==1
                index=index(idx);
                diameters(index,i)=data_current(k,1);
                data_new(index,:)=data_current(k,:);
            elseif length(FoundCrater)>1
                disp('length(FoundCrater)>1');
            end
        end
        if j>1
            m_data_pre=length(data_pre(:,1));
            for n=length(data_new(:,1)):-1:m_data_pre+1
                for k=n-1:-1:m_data_pre+1
                    if data_new(n,1)==data_new(k,1)&&data_new(n,2)==data_new(k,2)&&data_new(n,3)==data_new(k,3)
                        data_new(n,:)=[];
                        diameters(n,:)=[];
                        break;
                    end
                end
            end
        end
    end
    
    if ~isempty(data_pre)
        if (length(data_new(:,1))>length(data_pre(:,1)))
            data_pre(length(data_pre(:,1))+1:length(data_new(:,1)),:)=0;
        elseif (length(data_new(:,1))<length(data_pre(:,1)))
            data_new(length(data_new(:,1))+1:length(data_pre(:,1)),:)=0;
        end
    end
    data_pre=[data_pre,data_new];
    
    
end
for i=2:N
    index=diameters(:,1)>0&diameters(:,i)>0;
    systematic_error(i-1)=mean(diameters(index,i)./diameters(index,1))-1;
    
    if i==2
        idx=diameters(:,i)>0;
    else
        idx=idx&diameters(:,i)>0;
    end
end
systematic_error_diameter_mean=mean(systematic_error);

if N>2
    precision_diameter=std(diameters(idx,:),1,2);
else
    precision_diameter='N/A';
end
