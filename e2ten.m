function str=e2ten(str)
if nargin<1
    str='1.2e-05';
end
index_temp=strfind(str,'e');
index_temp=index_temp(end);
for j=index_temp+1:length(str)
    if strcmp(str(j),'-')
        continue;
    elseif str(j)>'0'
        break;
    else
        str(j)=[];
        break;
    end
end
str(index_temp:end+4)=['¡Á10^{',str(index_temp+1:end)];
str(end+1)='}';