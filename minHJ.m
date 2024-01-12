function [x,minf] = minHJ(f,a,b,eps)
format long;
if nargin < 1
    b1=3;
    b2=4;
    a=0.0027;
    ac=0.0026;
    a_eq=0.1;
    b_eq=1.83;
    f=@(D_mean)abs(log((a<=ac).*(a*(0.71*D_mean.^(1-0.036)).^-b2)+(a>ac).*(ac*(0.71*D_mean.^(1-0.036)).^-b2+(a-ac)*(0.71*D_mean.^(1-0.036)).^-b1))-log(a_eq*D_mean.^-b_eq));
    eps = 1.0e-14;
    a=0;
    b=100;
end
if nargin == 3
    eps = 1.0e-6;
end
l = a + 0.382*(b-a);
u = a + 0.618*(b-a);
k=1;
tol = b-a;
fl = f(l);
fu = f(u);
while tol>eps && k<100000
%     fl = subs(f , findsym(f), l);
%     fu = subs(f , findsym(f), u);
    
    if fl > fu
        a = l;
        l = u;
        u = a + 0.618*(b - a);
        
        fl = fu;
        fu = f(u);
    else
        b = u;
        u = l;
        l = a + 0.382*(b-a);
        
        fu = fl;
        fl = f(l);
    end
    k = k+1;
    tol = abs(b - a);
end
if k == 100000
    disp('找不到最小值！');
    x = NaN;
    minf = NaN;
    return;
end
x = (a+b)/2;
minf = f(x);
format short;

