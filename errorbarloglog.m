function hh = errorbarloglog(varargin)
%ERRORBAR Error bar plot.
%   ERRORBAR(X,Y,L,U) plots the graph of vector X vs. vector Y with
%   error bars specified by the vectors L and U.  L and U contain the
%   lower and upper error ranges for each point in Y.  Each error bar
%   is L(i) + U(i) long and is drawn a distance of U(i) above and L(i)
%   below the points in (X,Y).  The vectors X,Y,L and U must all be
%   the same length.  If X,Y,L and U are matrices then each column
%   produces a separate line.
%
%   ERRORBAR(X,Y,E) or ERRORBAR(Y,E) plots Y with error bars [Y-E Y+E].
%   ERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'.  The color is applied to the data line and
%   error bars while the linestyle and marker are applied to the data
%   line only.  See PLOT for possibilities.
%
%   ERRORBAR(AX,...) plots into AX instead of GCA.
%
%   H = ERRORBAR(...) returns a vector of errorbarseries handles in H.
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbar(x,y,e)
%   draws symmetric error bars of unit standard deviation.

%   L. Shure 5-17-88, 10-1-91 B.A. Jones 4-5-93
%   Copyright 1984-2015 MathWorks, Inc.

%---------------------------------------loglog (added by Xie)---------------%
% temp1=mfilename('fullpath');
% temp2=strfind(temp1,'\');
% path1=[temp1(1:temp2(end)) 'errorbar'];
% addpath(path1);
addpath('errorbar');
%-------------------------------------------------------------------------%

[~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);
nargs = length(args);
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end
[pvpairs,args,nargs,msg] = parseargs(args);
if ~isempty(msg), error(msg); end
if nargs < 2
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 4
    error(message('MATLAB:narginchk:tooManyInputs'));
end

hasXData = nargs ~= 2;
x = [];
switch nargs
    case 2
        [y,u] = deal(args{1:nargs});
        u = abs(u);
        l = u;
    case 3
        [x,y,u] = deal(args{1:nargs});
        if min(size(x))==1, x = x(:); end
        u = abs(u);
        l = u;
    case 4
        [x,y,l,u] = deal(args{1:nargs});
        if min(size(x))==1, x = x(:); end
end
if min(size(u))==1, u = u(:); end
if min(size(l))==1, l = l(:); end
if min(size(y))==1, y = y(:); end
n = size(y,2);
%---------------------------------------loglog (added by Xie)---------------%
temp=y<=l;
if sum(temp)>0
%     disp('Warning: lower error ranges for some points in Y larger than their corresponding Y,lower error ranges are assigned to be smaller than Y')
%     l(temp)=y(temp)-y(temp)*1e-15;
end
if sum(y<=0)
    disp('Warning: points in Y <= 0 are neglected')
end
%-------------------------------------------------------------------------%
% Make sure that x,y,l and u all are the same size:
if isempty(x)
    x = ones(size(y));
end
if ~isequal(size(x),size(y),size(l),size(u))
    error(message('MATLAB:errorbar:InputSizeMisMatch'));
end

% handle vectorized data sources and display names
extrapairs = cell(n,0);
if ~isempty(pvpairs) && (n > 1)
    [extrapairs, pvpairs] = vectorizepvpairs(pvpairs,n,...
        {'XDataSource','YDataSource',...
        'UDataSource','LDataSource',...
        'DisplayName'});
end

if isempty(cax) || ishghandle(cax,'axes')
%     cax = newplot(cax);%commended by Xie---------------%
    cax = gca;%commended by Xie---------------%
    parax = cax;
    hold_state = ishold(cax);
else
    parax = cax;
    cax = ancestor(cax,'axes');
    hold_state = true;
end

h = [];
colorPropName = 'Color';
autoColor = ~any(strcmpi('color',pvpairs(1:2:end)));
if autoColor
    colorPropName = 'Color_I';
end
stylePropName = 'LineStyle';
autoStyle = ~any(strcmpi('linestyle',pvpairs(1:2:end)));
if autoStyle
    stylePropName = 'LineStyle_I';
end
xdata = {};
for k=1:n
    % extract data from vectorizing over columns
    if hasXData
        xdata = {'XData', datachk(x(:,k))};
    end
    [ls,c,m] = nextstyle(cax,autoColor,autoStyle);
    if k==1
        h = matlab.graphics.chart.primitive.ErrorBar('YData',datachk(y(:,k)),...
            'UData',datachk(u(:,k)),...
            'LData',datachk(l(:,k)),xdata{:},...
            colorPropName,c,stylePropName,ls,'Marker_I',m,...
            pvpairs{:},extrapairs{k,:},'Parent',parax);
    else
        h(k) = matlab.graphics.chart.primitive.ErrorBar('YData',datachk(y(:,k)),...
            'UData',datachk(u(:,k)),...
            'LData',datachk(l(:,k)),xdata{:},...
            colorPropName,c,stylePropName,ls,'Marker_I',m,...
            pvpairs{:},extrapairs{k,:},'Parent',parax);
    end
end
if ~hold_state
    set(cax,'Box','on');
end
if nargout>0, hh = h; end
%---------------------------------------loglog (added by Xie)---------------%
set(gca,'yscale','log')
set(gca,'xscale','log')
% rmpath(path1);
rmpath('errorbar');
%-------------------------------------------------------------------------%
function [pvpairs,args,nargs,msg] = parseargs(args)
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);
% check for LINESPEC
if ~isempty(pvpairs)
    [l,c,m,tmsg]=colstyle(pvpairs{1},'plot');
    if isempty(tmsg)
        pvpairs = pvpairs(2:end);
        if ~isempty(l)
            pvpairs = [{'LineStyle',l},pvpairs];
        end
        if ~isempty(c)
            pvpairs = [{'Color',c},pvpairs];
        end
        if ~isempty(m)
            pvpairs = [{'Marker',m},pvpairs];
        end
    end
end
msg = checkpvpairs(pvpairs);
nargs = length(args);
