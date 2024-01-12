function str=formatfloat(x)
if nargin<1
    x=-10;
end
e=log10(abs(x));
e=floor(e);

c=x./10.^e;
for i=1:length(x)
    if e~=0
        str=sprintf('%.1f¡Á10^{%d}',c(i),e(i));
    else
        str=sprintf('%.1f',c(i));
    end
end