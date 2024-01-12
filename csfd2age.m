function [age,DL_fit,cumN_fit]=csfd2age(DL,cumN,err,area,PF_type)
%[age,DL_fit,cumN_fit]=csfd2age(DL,cumN,err,area)
if nargin<1
    close all
    folder='I:\Documents\Paper\Small craters\ArcGIS\Serenitatis\»®·ÖµØ²ãV1\shp\Area_i\';
    filename_area=[folder  '4.dbf'];
    data=shaperead(filename_area);
    area=data.Area;
    filename=[folder 'Craters_4.dbf'];
    [DL,cumN,err]=shp2cum(filename,area,'pseudo-log');%2^(1/8)
    errorbarloglog(DL,cumN,err,'ro');hold on
    
    PF_type='NPF';
     PF_type='Xie_2017';
    PF_type='b=3';%N=D^-b
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
[N1_fitting,N1_err,DL_fit,cumN_fit]=csfd2density(DL,cumN,err,area,PF_type);
age=density2age(N1_fitting);

age_lower=density2age(N1_fitting-N1_err);
age_upper=density2age(N1_fitting+N1_err);

age(2)=age_lower;
age(3)=age_upper;

if nargin<1
    loglog(DL_fit,cumN_fit,'r-');
end



