function [EPF_dsfd,EPF_isfd,EPF_csfd]=diameter2edf(D,alpha,edf_D)
if nargin<1
    D=1+rand(1000,1)*100;
end
if nargin<2
    alpha=0.1;
end
if nargin<3
    edf_D=logspace(log10(min(D)/2),log10(max(D)*2),1000);
end
EPF_dsfd=zeros(size(edf_D));
step=edf_D(2)/edf_D(1);
sigma=D*alpha;
for i=1:length(D)
    index=find(D(i)-4*sigma(i)<edf_D*step&edf_D<=D(i)+4*sigma(i));
    y=1/(sqrt(2*pi)*sigma(i))*exp(-0.5*((edf_D(index)-D(i))/sigma(i)).^2);
%     y=y/sum(y.*bin_width(index))*MOPF_isfd(i);%%不能通过除以sum(y.*bin_width(index))来归一化，因为有可能某些情况是在横坐标范围外的，比如L=0.1,而一个0.1直径的撞击坑，有一半概率是小于L的
    EPF_dsfd(index)=EPF_dsfd(index)+y;
%     p(i,index)=y;
%     loglog(edf.DL(index),y,'r--');hold on
end

if nargout>1||nargin<1
    edf_DR=edf_D*step^0.5;
    edf_DL=edf_D*step^-0.5;
    bin_width=edf_DR-edf_DL;
    EPF_isfd=EPF_dsfd.*bin_width;
end

if nargout>2||nargin<1
    EPF_csfd=cumsum(EPF_isfd(end:-1:1));
    EPF_csfd(end:-1:1)=EPF_csfd;
end

if nargin<1
    loglog(edf_D,EPF_dsfd,'r',edf_D,EPF_isfd,'g',edf_D,EPF_csfd,'b');
    legend('DSFD','ISFD','CSFD')
end